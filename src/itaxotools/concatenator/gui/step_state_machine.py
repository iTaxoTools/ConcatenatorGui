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

from . import step_progress_bar


class NavigateEvent(QtCore.QEvent):
    """Custom event for use in StepStateMachines"""
    userEvent = QtCore.QEvent.registerEventType()
    events = set()

    class Action(enum.Enum):
        Back = enum.auto()
        Next = enum.auto()
        Exit = enum.auto()
        Cancel = enum.auto()
        New = enum.auto()

    def __init__(self, where: Action):
        """Pass name and args"""
        super().__init__(QtCore.QEvent.Type(self.userEvent))
        self.where = where
        # Avoid garbage-collection
        NavigateEvent.events.add(self)


class NavigateTransition(QtStateMachine.QAbstractTransition):
    """Custom transition for use in StepStateMachines"""

    def __init__(self, where: NavigateEvent.Action):
        """Only catch events with given direction"""
        super().__init__()
        self.where = where

    def eventTest(self, event):
        """Check for NavigateEvent"""
        if event.type() == NavigateEvent.userEvent:
            return event.where == self.where
        return False

    def onTransition(self, event):
        """Override virtual function"""
        # Allow event to be garbage-collected
        # NamedEvent.events.remove(event)


class StepState(QtStateMachine.QState):
    def __init__(self, name, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._nextState = None
        self._prevState = None
        self.stepProgressBar = self.parent().stepProgressBar
        self.navigationFooter = self.parent().navigationFooter
        self.name = name

        self._nextTransition = NavigateTransition(NavigateEvent.Action.Next)
        self.addTransition(self._nextTransition)
        self._prevTransition = NavigateTransition(NavigateEvent.Action.Back)
        self.addTransition(self._prevTransition)

    def onEntry(self, event):
        super().onEntry(event)
        self.stepProgressBar.activateKey(self.name)
        if self.nextState:
            if self.prevState:
                self.navigationFooter.showIntermediate()
            else:
                self.navigationFooter.showBegin()
        else:
            self.navigationFooter.showFinal()

    @property
    def nextState(self):
        return self._nextState

    @nextState.setter
    def nextState(self, state):
        self._nextTransition.setTargetState(state)
        self._nextState = state

    @property
    def prevState(self):
        return self._prevState

    @prevState.setter
    def prevState(self, state):
        self._prevTransition.setTargetState(state)
        self._prevState = state


class StepStateMachine(QtStateMachine.QStateMachine):
    """Provides convenience functionality for StepStates"""

    def __init__(self,
                 parent: QtWidgets.QWidget,
                 stepProgressBar: step_progress_bar.StepProgressBar,
                 navigationFooter: widgets.NavigationFooter
                 ):
        super().__init__(parent)
        self.stepProgressBar = stepProgressBar
        self.navigationFooter = navigationFooter
        self.states = {}
        self.steps = []

        self.navigationFooter.setButtonActions({
            'next': self.eventGenerator(action=NavigateEvent.Action.Next),
            'back': self.eventGenerator(action=NavigateEvent.Action.Back),
            'exit': self.eventGenerator(action=NavigateEvent.Action.Exit),
            'cancel': self.eventGenerator(action=NavigateEvent.Action.Cancel),
            'new': self.eventGenerator(action=NavigateEvent.Action.New),
        })

    def eventGenerator(self, checked=False, action=None):
        def eventFunction():
            self.postEvent(NavigateEvent(action))
        return eventFunction

    def addStep(self, name, text, weight=1, visible=True):
        state = StepState(name, self)

        if len(self.steps) == 0:
            self.setInitialState(state)
            state.prevState = None
        else:
            prev = self.steps[-1]
            state.prevState = prev
            prev.nextState = state

        self.stepProgressBar.addStep(name, text, weight, visible)
        self.states[name] = state
        self.steps.append(state)
