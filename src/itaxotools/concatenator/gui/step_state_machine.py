# -----------------------------------------------------------------------------
# ConcatenatorQt - GUI for Concatenator
# Copyright (C) 2021  Patmanidis Stefanos
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------

"""State machine logic for a series of steps"""


from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtStateMachine

import enum

from itaxotools.common import widgets

from . import step_progress_bar as spb


class FailedError(Exception):
    def __init__(self, code=1):
        super().__init__(f'Thread failed with code: {code}')
        self.code = code


class CancelledError(Exception):
    def __init__(self, code=-1):
        super().__init__(f'Thread cancelled with code: {code}')
        self.code = code


class TerminatedError(Exception):
    def __init__(self):
        super().__init__('Thread was forcibly terminated')


class WorkerThread(QtCore.QThread):
    done = QtCore.Signal(object)
    fail = QtCore.Signal(object)
    cancel = QtCore.Signal(object)

    def __init__(self, function, *args, **kwargs):
        super().__init__()
        self.function = function
        self.args = args
        self.kwargs = kwargs

    def run(self):
        try:
            result = self.function(*self.args, **self.kwargs)
        except CancelledError as exception:
            self.cancel.emit(exception)
        except Exception as exception:
            self.fail.emit(exception)
        else:
            self.done.emit(result)


class NavigateEvent(QtCore.QEvent):
    """Custom event for use in StepStateMachines"""
    userEvent = QtCore.QEvent.registerEventType()
    events = set()

    class Event(enum.Enum):
        Back = enum.auto()
        Next = enum.auto()
        Exit = enum.auto()
        Done = enum.auto()
        Fail = enum.auto()
        Cancel = enum.auto()
        New = enum.auto()

    def __init__(self, event: Event):
        """Pass name and args"""
        super().__init__(QtCore.QEvent.Type(self.userEvent))
        self.event = event
        # ! Avoid garbage-collection
        NavigateEvent.events.add(self)


class NavigateTransition(QtStateMachine.QAbstractTransition):
    """Custom transition for use in StepStateMachines"""

    def __init__(self, event: NavigateEvent.Event, filter=lambda e: True):
        """Only catch events with given direction"""
        super().__init__()
        self.event = event
        self.filter = filter

    def eventTest(self, event):
        """Check for NavigateEvent"""
        if event.type() == NavigateEvent.userEvent:
            if event.event == self.event:
                return self.filter(event)
        return False

    def onTransition(self, event):
        """Override virtual function"""
        # ! Allow event to be garbage-collected
        QtCore.QTimer.singleShot(0, lambda: NavigateEvent.events.remove(event))


class StepSubState(QtStateMachine.QState):
    def __init__(self, parent, stack, navMode, progStatus):
        super().__init__(parent)
        self.stack = stack
        self.navMode = navMode
        self.progStatus = progStatus
        self.stepProgressBar = parent.stepProgressBar
        self.navigationFooter = parent.navigationFooter
        self.widget = self.draw()
        self.stack.addWidget(self.widget)

    def draw(self):
        return QtWidgets.QWidget()

    def onEntry(self, event):
        backwards = False
        if isinstance(event, NavigateEvent):
            backwards = event.event == NavigateEvent.Event.Back
        self.stack.setCurrentWidget(self.widget)
        self.navigationFooter.setMode(self.navMode, backwards)
        self.stepProgressBar.setStatus(self.progStatus)


class StepState(StepSubState):
    def __init__(self, name, parent, stack):
        super().__init__(parent, stack, None, spb.states.Active)
        self.nextState = None
        self.prevState = None
        self.name = name
        self.cog()

    def cog(self):
        transitions = {}
        transitions['next'] = NavigateTransition(NavigateEvent.Event.Next)
        transitions['back'] = NavigateTransition(NavigateEvent.Event.Back)
        for transition in transitions:
            self.addTransition(transitions[transition])
        self.transitions = transitions

    def setNextState(self, state):
        self.transitions['next'].setTargetState(state)
        self.nextState = state

    def setPrevState(self, state):
        self.transitions['back'].setTargetState(state)
        self.prevState = state

    def onEntry(self, event):
        if self.nextState:
            if self.prevState:
                self.navMode = widgets.NavigationFooter.Mode.Middle
            else:
                self.navMode = widgets.NavigationFooter.Mode.First
        else:
            self.navMode = widgets.NavigationFooter.Mode.Final
        self.stepProgressBar.activateKey(self.name)
        super().onEntry(event)


class StepTriState(StepState):

    StepEdit = StepSubState
    StepWait = StepSubState
    StepFail = StepSubState

    class DataObject(object):
        pass

    def __init__(self, name, parent, stack):
        super().__init__(name, parent, stack)
        parent.signalTerminate.connect(self.threadTerminate)
        self.workerThread = WorkerThread(self.work)
        self.workerThread.done.connect(self._onDone)
        self.workerThread.fail.connect(self._onFail)
        self.workerThread.cancel.connect(self.onCancel)
        self.threadTerminateTimeout = 1000
        self.data = self.DataObject()

    def draw(self):
        widget = QtWidgets.QStackedWidget()
        widget.setContentsMargins(0, 0, 0, 0)
        return widget

    def threadTerminate(self):
        """Attempt to cleanly exit, then forcibly terminate (dangerous!)"""
        self.parent().waiting = False
        self.workerThread.requestInterruption()
        self.workerThread.exit(-1)
        if not self.workerThread.wait(self.threadTerminateTimeout):
            self.workerThread.terminate()
            self.workerThread.wait()
            self.onCancel(TerminatedError())

    def threadExit(self, code: int = 0):
        """
        Call this from within a worker thread to indicate the work is done.
        Argument `code` determines transition behaviour: zero for success,
        a positive value for failure, a negative value on interruption.
        """
        self.workerThread.exit(code)

    def work(self):
        """
        This is executed with StepWait.onEntry on a new QThread.
        Overload to assign a meaningful task. This should either:
        - Return a result (forwarded to self.onDone()),
        - Raise an exception (forwarded to self.onFail() or self.onCancel()),
        - Queue self.threadExit() and call super().work().
        For long tasks, you may check workerThread.isInterruptionRequested()
        and raise CancelledError() if True.
        """
        code = self.thread().exec()
        if code < 0:
            raise CancelledError(code)
        if code > 0:
            raise FailedError(code)

    def _onDone(self, result):
        self.onDone(result)
        self.machine().postEvent(NavigateEvent(NavigateEvent.Event.Done))

    def _onFail(self, exception):
        self.onFail(exception)
        self.machine().postEvent(NavigateEvent(NavigateEvent.Event.Fail))

    def onDone(self, result):
        """Called after work() is successfully completed"""
        pass

    def onFail(self, exception):
        """Called after work() has raised an exception"""
        pass

    def onCancel(self, exception):
        """Called after work() was cancelled"""
        pass

    def filterNext(self):
        """Intercepts transition editNext"""
        return True

    def filterBack(self):
        """Intercepts transition dataBack"""
        return True

    def filterCancel(self):
        """Intercepts transition waitCancel"""
        return True

    def wrapStepWait(self):
        """Wrap self.StepWait methods to handle worker"""
        tristate = self

        class _StepWait(self.StepWait):
            def onEntry(self, event):
                tristate.parent().waiting = True
                tristate.workerThread.start()
                super().onEntry(event)

            def onExit(self, event):
                tristate.parent().waiting = False
                tristate.threadTerminate()
                super().onExit(event)
        return _StepWait

    def cog(self):
        self.states = {}
        nm = widgets.NavigationFooter.Mode
        ps = spb.states
        self.addSubState('edit', self.StepEdit, nm.Middle, ps.Active)
        self.addSubState('fail', self.StepFail, nm.Error, ps.Failed)
        self.addSubState('wait', self.wrapStepWait(), nm.Wait, ps.Ongoing)
        self.setInitialState(self.states['edit'])

        ev = NavigateEvent.Event
        self.transitions = {}
        self.addSubTransition(
            'editNext', ev.Next, lambda e: self.filterNext(),
            self.states['edit'], self.states['wait'])
        self.addSubTransition(
            'waitFail', ev.Fail, lambda e: True,
            self.states['wait'], self.states['fail'])
        self.addSubTransition(
            'waitCancel', ev.Cancel, lambda e: self.filterCancel(),
            self.states['wait'], self.states['edit'])
        self.addSubTransition(
            'failBack', ev.Back, lambda e: True,
            self.states['fail'], self.states['edit'])
        self.addSubTransition(
            'waitDone', ev.Done, lambda e: True,
            self.states['wait'], None)
        self.addSubTransition(
            'dataBack', ev.Back, lambda e: self.filterBack(),
            self.states['edit'], None)

    def addSubState(self, name, cls, navMode, progStatus):
        self.states[name] = cls(self, self.widget, navMode, progStatus)

    def addSubTransition(self, name, action, filter, source, target):
        transition = NavigateTransition(action, filter)
        transition.setTargetState(target)
        source.addTransition(transition)
        self.transitions[name] = transition

    def setNextState(self, state):
        self.transitions['waitDone'].setTargetState(state)
        self.nextState = state

    def setPrevState(self, state):
        self.transitions['dataBack'].setTargetState(state)
        self.prevState = state


class StepStateMachine(QtStateMachine.QStateMachine):
    """Provides convenience functionality for StepStates"""

    signalTerminate = QtCore.Signal()

    def __init__(self,
                 parent: QtWidgets.QWidget,
                 stepProgressBar: spb.StepProgressBar,
                 navigationFooter: widgets.NavigationFooter,
                 stack: QtWidgets.QStackedLayout
                 ):
        super().__init__(parent)
        self.stepProgressBar = stepProgressBar
        self.navigationFooter = navigationFooter
        self.stack = stack
        self.waiting = False
        self.states = {}
        self.steps = []

        self.navigationFooter.setButtonActions({
            'next': self.eventGenerator(action=NavigateEvent.Event.Next),
            'back': self.eventGenerator(action=NavigateEvent.Event.Back),
            'exit': self.eventGenerator(action=NavigateEvent.Event.Exit),
            'cancel': self.eventGenerator(action=NavigateEvent.Event.Cancel),
            'new': self.eventGenerator(action=NavigateEvent.Event.New),
        })

    def eventGenerator(self, checked=False, action=None):
        # This will be called by QPushButton.clicked,
        # so `checked` will be the first kwarg.
        def eventFunction():
            self.postEvent(NavigateEvent(action))
        return eventFunction

    def addStep(self, name, text=None, weight=1, visible=True, cls=StepState):
        text = name if text is None else text
        state = cls(name, self, self.stack)

        if len(self.steps) == 0:
            self.setInitialState(state)
            state.prevState = None
        else:
            prev = self.steps[-1]
            state.setPrevState(prev)
            prev.setNextState(state)

        self.stepProgressBar.addStep(name, text, weight, visible)
        self.states[name] = state
        self.steps.append(state)

    def terminate(self):
        """Attempt to cleanly exit child-states, then forcibly terminate"""
        self.signalTerminate.emit()

    def cancel(self):
        """If in wait check, try to cancel, else return None"""
        if self.waiting:
            self.postEvent(NavigateEvent(NavigateEvent.Event.Cancel))
            return True
        else:
            return None
