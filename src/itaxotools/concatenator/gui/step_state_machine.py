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


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__dict__ = self


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
        self.timeout = 1000

    def run(self):
        try:
            result = self.function(*self.args, **self.kwargs)
        except CancelledError as exception:
            self.cancel.emit(exception)
        except Exception as exception:
            self.fail.emit(exception)
        else:
            self.done.emit(result)

    def check(self):
        """
        Call this from within a worker thread to check if the thread
        should exit by user request. If so, CancelledError is raised,
        which should then be handled by WorkerThread.run().
        """
        if self.isInterruptionRequested():
            raise CancelledError(-1)

    def exit(self, code: int = 0):
        """
        Call this from within a worker thread to indicate the work is done.
        Argument `code` determines transition behaviour: zero for success,
        a positive value for failure, a negative value on interruption.
        """
        super().exit(code)

    def terminate(self):
        """
        Attempt to cleanly exit, then forcibly terminate
        after a short timeout (dangerous!)
        """
        self.requestInterruption()
        self.exit(-1)
        if not self.wait(self.timeout):
            super().terminate()
            self.wait()
            self.cancel.emit(TerminatedError())


class NavigateAction(enum.Enum):
    Back = enum.auto()
    Next = enum.auto()
    Exit = enum.auto()
    Done = enum.auto()
    Fail = enum.auto()
    Cancel = enum.auto()
    New = enum.auto()
    Skip = enum.auto()


class NavigateEvent(QtStateMachine.QStateMachine.SignalEvent):
    """
    Wraps a SignalEvent meant for a NavigateTransition.
    The signal must be emitted with a NavigateAction as the first argument.
    """

    def __init__(self, event):
        if self.isValidCast(event):
            self._action = event.arguments()[0]
            arguments = event.arguments()[1:]
            super().__init__(event.sender(), event.signalIndex(), arguments)
        else:
            self._action = None
            super().__init__(None, -1, [])
            self.type = event.type

    def action(self):
        """Get the associated NavigateAction"""
        return self._action

    @staticmethod
    def isValidCast(event):
        return (isinstance(event, QtStateMachine.QStateMachine.SignalEvent)
                and event.arguments()
                and isinstance(event.arguments()[0], NavigateAction))


class NavigateTransition(QtStateMachine.QSignalTransition):
    """Only catches NavigateEvents of a specific NavigateAction."""

    def __init__(
        self, signal,
        action: NavigateAction,
        filter=lambda e: True,
        function=lambda: None,
    ):
        super().__init__(signal)
        self.action = action
        self.filter = filter
        self.function = function

    def eventTest(self, event):
        action = NavigateEvent(event).action()
        if action == self.action:
            return self.filter(event)
        return False

    def onTransition(self, event):
        self.function()
        super().onTransition(event)


class StepSubState(QtStateMachine.QState):

    title = None
    description = None

    def __init__(self, parent, stack):
        super().__init__(parent)
        self.stack = stack
        self.navMode: widgets.NavigationFooter.Mode = None
        self.progStatus: spb.states.AbstractStatus = None
        self.stepProgressBar = parent.stepProgressBar
        self.header = parent.header
        self.footer = parent.footer
        self.widget = self.draw()
        self.stack.addWidget(self.widget)
        self.header.showTask(self.title, self.description)
        self.header.updateLabelWidth()

    def draw(self):
        return QtWidgets.QWidget()

    def onEntry(self, event):
        backwards = False
        action = NavigateEvent(event).action()
        backwards = action == NavigateAction.Back
        self.stack.setCurrentWidget(self.widget)
        if self.title is None and self.description is None:
            self.header.showTool()
        else:
            self.header.showTask(self.title, self.description)
        if self.navMode:
            self.footer.setMode(self.navMode, backwards)
        if self.progStatus:
            self.stepProgressBar.setStatus(self.progStatus)

    def clear(self):
        """Virtual. Override should reset local data and widgets"""
        pass


class StepState(StepSubState):
    def __init__(self, name, parent, stack):
        super().__init__(parent, stack)
        self.progStatus = spb.states.Active
        self.nextState = None
        self.prevState = None
        self.name = name
        self.cog()

    def cog(self):
        m = self.machine()
        transitions = AttrDict()
        transitions.next = m.navigateTransition(NavigateAction.Next)
        transitions.back = m.navigateTransition(NavigateAction.Back)
        transitions.exit = m.navigateTransition(NavigateAction.Exit)
        transitions.exit.setTargetState(m.states.final)
        for transition in transitions:
            self.addTransition(transitions[transition])
        self.transitions = transitions

    def setNextState(self, state):
        self.transitions.next.setTargetState(state)
        self.transitions.exit.setTargetState(None)
        self.nextState = state

    def setPrevState(self, state):
        self.transitions.back.setTargetState(state)
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


class StepTriStateEdit(StepSubState):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.navMode = widgets.NavigationFooter.Mode.Middle
        self.progStatus = spb.states.Active


class StepTriStateWait(StepSubState):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.navMode = widgets.NavigationFooter.Mode.Wait
        self.progStatus = spb.states.Ongoing
        self.parent().updated.connect(self._update)

    def onEntry(self, event):
        self.parent().worker.start()
        super().onEntry(event)

    def onExit(self, event):
        self.parent().worker.terminate()
        super().onExit(event)

    def _update(self, args, kwargs):
        self.update(*args, **kwargs)

    def update(self, *args, **kwargs):
        """Virtual. Called when worker thread uses update()"""
        pass


class StepTriStateDone(StepSubState):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.navMode = widgets.NavigationFooter.Mode.Middle
        self.progStatus = spb.states.Complete
        self.result = None

    def draw(self):
        return self.parent().states.wait.widget


class StepTriStateFail(StepSubState):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.navMode = widgets.NavigationFooter.Mode.Error
        self.progStatus = spb.states.Failed
        self.exception = None

    def draw(self):
        return self.parent().states.wait.widget


class StepTriState(StepState):

    StepEdit = StepTriStateEdit
    StepWait = StepTriStateWait
    StepDone = StepTriStateDone
    StepFail = StepTriStateFail

    updated = QtCore.Signal(list, dict)

    class DataObject(object):
        pass

    def __init__(self, name, parent, stack):
        super().__init__(name, parent, stack)
        parent.installWorker(self)
        self.data = self.DataObject()

    def draw(self):
        widget = QtWidgets.QStackedWidget()
        widget.setContentsMargins(0, 0, 0, 0)
        return widget

    def work(self):
        """
        Virtual. Executed by StepWait.onEntry on a new QThread.
        Overload to assign a meaningful task. This should either:
        - Return a result (forwarded to self.onDone()),
        - Raise an exception (forwarded to self.onFail() or self.onCancel()),
        - Queue self.worker.exit() and call super().work().
        For long tasks, you may use self.worker.check() at regular
        intervals to cleanly interrupt execution at user request.
        You may use self.update() to update StepWait contents.
        """
        code = self.worker.exec()
        if code < 0:
            raise CancelledError(code)
        if code > 0:
            raise FailedError(code)

    def update(self, *args, **kwargs):
        self.updated.emit(args, kwargs)

    def onEntry(self, event):
        if self.skipAll():
            action = NavigateEvent(event).action()
            self.machine().navigate(action)
        else:
            super().onEntry(event)

    def _onDone(self, result):
        self.onDone(result)
        self.states.done.result = result
        self.states.fail.exception = None
        self.machine().navigate(NavigateAction.Done)

    def _onFail(self, exception):
        self.onFail(exception)
        self.states.done.result = None
        self.states.fail.exception = exception
        self.machine().navigate(NavigateAction.Fail)

    def onDone(self, result):
        """Called after work() is successfully completed"""
        pass

    def onFail(self, exception):
        """Called after work() has raised an exception"""
        pass

    def onCancel(self, exception):
        """Called after work() was cancelled"""
        pass

    def _skip(self):
        self.machine().navigate(NavigateAction.Skip)

    def skipAll(self):
        """Override to skip the whole step"""
        return False

    def skipWait(self):
        """Override to skip the waiting step"""
        return False

    def _filterNext(self, event):
        if self.skipAll():
            self._skip()
            return False
        if self.skipWait():
            if self.filterSkip(event):
                self._skip()
            return False
        return self.filterNext(event)

    def filterNext(self, event):
        """Intercepts transition editNext"""
        return True

    def filterSkip(self, event):
        """Intercepts skipWait effects"""
        return True

    def filterBack(self, event):
        """Intercepts transition dataBack"""
        return True

    def filterCancel(self, event):
        """Intercepts transition waitCancel"""
        return True

    def cog(self):
        self.states = AttrDict()
        self.addSubState('edit', self.StepEdit)
        self.addSubState('wait', self.StepWait)
        self.addSubState('done', self.StepDone)
        self.addSubState('fail', self.StepFail)
        self.setInitialState(self.states.edit)

        self.transitions = AttrDict()
        self.addSubTransition(
            'editNext', NavigateAction.Next,
            self.states.edit, self.states.wait,
            self._filterNext)
        self.addSubTransition(
            'editSkip', NavigateAction.Skip,
            self.states.edit, None)
        self.addSubTransition(
            'waitFail', NavigateAction.Fail,
            self.states.wait, self.states.fail)
        self.addSubTransition(
            'waitDone', NavigateAction.Done,
            self.states.wait, self.states.done)
        self.addSubTransition(
            'waitCancel', NavigateAction.Cancel,
            self.states.wait, self.states.edit,
            self.filterCancel)
        self.addSubTransition(
            'failBack', NavigateAction.Back,
            self.states.fail, self.states.edit)
        self.addSubTransition(
            'doneBack', NavigateAction.Back,
            self.states.done, self.states.edit)
        self.addSubTransition(
            'doneNext', NavigateAction.Next,
            self.states.done, None)
        self.addSubTransition(
            'editBack', NavigateAction.Back,
            self.states.edit, None,
            self.filterBack)

    def addSubState(self, name, cls):
        self.states[name] = cls(self, self.widget)

    def addSubTransition(
        self, name, action, source, target,
        filter=lambda e: True
    ):
        transition = self.machine().navigateTransition(action, filter)
        transition.setTargetState(target)
        source.addTransition(transition)
        self.transitions[name] = transition

    def setNextState(self, state):
        self.transitions.doneNext.setTargetState(state)
        self.transitions.editSkip.setTargetState(state)
        self.nextState = state

    def setPrevState(self, state):
        self.transitions.editBack.setTargetState(state)
        self.prevState = state

    def clear(self):
        for state in dict(self.states).values():
            state.clear()


class StepStateMachine(QtStateMachine.QStateMachine):
    """Provides extra functionality for StepStates"""

    terminateSignal = QtCore.Signal()
    navigateSignal = QtCore.Signal(NavigateAction)

    def __init__(
        self,
        parent: QtWidgets.QWidget,
        stepProgressBar: spb.StepProgressBar,
        header: widgets.Header,
        footer: widgets.NavigationFooter,
        stack: QtWidgets.QStackedLayout
     ):
        super().__init__(parent)
        self.stepProgressBar = stepProgressBar
        self.header = header
        self.footer = footer
        self.stack = stack
        self.states = AttrDict()
        self.steps = []
        self.waiting = False

        self.states.final = QtStateMachine.QFinalState(self)

        self.footer.setButtonActions({
            'next': self.eventGenerator(action=NavigateAction.Next),
            'back': self.eventGenerator(action=NavigateAction.Back),
            'exit': self.eventGenerator(action=NavigateAction.Exit),
            'cancel': self.eventGenerator(action=NavigateAction.Cancel),
            'new': self.eventGenerator(action=NavigateAction.New),
        })

    @QtCore.Slot()
    def _setWaiting(self):
        self.waiting = True

    @QtCore.Slot()
    def _unsetWaiting(self):
        self.waiting = False

    def installWorker(self, state):
        """Add a WorkerThread to target state and configure signals"""
        state.worker = WorkerThread(state.work)
        onDone = getattr(state, '_onDone', state.onDone)
        onFail = getattr(state, '_onFail', state.onFail)
        onCancel = getattr(state, '_onCancel', state.onCancel)
        state.worker.done.connect(onDone)
        state.worker.fail.connect(onFail)
        state.worker.cancel.connect(onCancel)
        state.worker.started.connect(self._setWaiting)
        state.worker.finished.connect(self._unsetWaiting)
        self.terminateSignal.connect(state.worker.terminate)

    def eventGenerator(self, checked=False, action=None):
        # This will be called by QPushButton.clicked,
        # so `checked` will be the first kwarg.
        def eventFunction():
            self.navigate(action)
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
        self.terminateSignal.emit()

    def navigate(self, action):
        self.navigateSignal.emit(action)

    def navigateTransition(
            self, action,
            filter=lambda e: True,
            function=lambda: None):
        return NavigateTransition(
            self.navigateSignal, action, filter, function)

    def navigateTransitionClear(self):
        return self.navigateTransition(
            NavigateAction.New, function=self.clear)

    def clear(self):
        """Clears local data for all states"""
        for state in dict(self.states).values():
            if isinstance(state, StepSubState):
                state.clear()

    def cancel(self):
        """If in wait check, try to cancel"""
        if self.waiting:
            self.navigate(NavigateAction.Cancel)
            return False
        else:
            return True
