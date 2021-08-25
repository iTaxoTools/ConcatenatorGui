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

    def __init__(self, event: NavigateEvent.Event):
        """Only catch events with given direction"""
        super().__init__()
        self.event = event

    def eventTest(self, event):
        """Check for NavigateEvent"""
        if event.type() == NavigateEvent.userEvent:
            return event.event == self.event
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
        super().__init__(parent, stack, None, None)
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
        self.progStatus = spb.states.Active
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

    # Overload these for custom sub-states
    StepEdit = StepSubState
    StepWait = StepSubState
    StepFail = StepSubState

    def draw(self):
        return QtWidgets.QStackedWidget()

    def cog(self):
        self.states = {}
        nm = widgets.NavigationFooter.Mode
        ps = spb.states
        self.addSubState('edit', self.StepEdit, nm.Middle, ps.Active)
        self.addSubState('wait', self.StepWait, nm.Wait, ps.Ongoing)
        self.addSubState('fail', self.StepFail, nm.Error, ps.Failed)
        self.setInitialState(self.states['edit'])

        ev = NavigateEvent.Event
        self.transitions = {}
        self.addSubTransition(
            'editNext', self.states['edit'], self.states['wait'], ev.Next)
        self.addSubTransition(
            'waitFail', self.states['wait'], self.states['fail'], ev.Fail)
        self.addSubTransition(
            'waitCancel', self.states['wait'], self.states['edit'], ev.Cancel)
        self.addSubTransition(
            'failBack', self.states['fail'], self.states['edit'], ev.Back)
        self.addSubTransition(
            'waitDone', self.states['wait'], None, ev.Done)
        self.addSubTransition(
            'dataBack', self.states['edit'], None, ev.Back)

    def addSubState(self, name, cls, navMode, progStatus):
        self.states[name] = cls(self, self.widget, navMode, progStatus)

    def addSubTransition(self, name, source, target, action):
        transition = NavigateTransition(action)
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
