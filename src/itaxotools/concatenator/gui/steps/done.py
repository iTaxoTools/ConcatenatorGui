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

"""StepDone"""

from PySide6 import QtCore
from PySide6 import QtWidgets

from .. import step_state_machine as ssm

from itaxotools import common
import itaxotools.common.resources # noqa


class DataObject(object):
    pass


class StepDone(ssm.StepState):

    title = 'Concatenation Complete'
    description = 'Successfully exported output'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = DataObject()

    def cog(self):
        super().cog()
        machine = self.machine()
        transition = machine.navigateTransitionClear()
        transition.setTargetState(machine.states.input)
        self.addTransition(transition)
        self.transitions.new = transition

    def draw(self):
        widget = QtWidgets.QWidget()

        path = common.resources.get(
            'itaxotools.concatenator.gui', 'docs/done.md')
        with open(path) as file:
            text = file.read()

        self.label = QtWidgets.QLabel()
        self.label.setTextFormat(QtCore.Qt.MarkdownText)
        self.label.setOpenExternalLinks(True)
        self.label.setText(text)
        self.label.setStyleSheet("background: transparent; border: 0px;")

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.label)
        layout.addStretch(1)
        layout.setSpacing(24)
        layout.setContentsMargins(0, 0, 0, 32)
        widget.setLayout(layout)

        return widget

    def onEntry(self, event):
        super().onEntry(event)
