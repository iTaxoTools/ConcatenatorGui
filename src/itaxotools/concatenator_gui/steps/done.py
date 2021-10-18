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

from itaxotools import common
import itaxotools.common.resources # noqa

from .. import step_state_machine as ssm
from .. import widgets

from . import wait


class StepDone(ssm.StepState):

    title = 'Results exported successfully'
    description = 'Concatenation complete'

    def cog(self):
        super().cog()
        machine = self.machine()
        transition = machine.navigateTransitionClear()
        transition.setTargetState(machine.states.input)
        self.addTransition(transition)
        self.transitions.new = transition

    def draw(self):
        widget = QtWidgets.QWidget()

        self.progress = wait.ProgressBar()
        self.progress.bar.setMaximum(1)
        self.progress.bar.setValue(1)

        self.confirm = self.progress.label
        self.confirm.setTextFormat(QtCore.Qt.RichText)
        self.confirm.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)

        path = common.resources.get(
            'itaxotools.concatenator_gui', 'docs/done.html')
        self.label = widgets.HtmlLabel(path)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.progress)
        layout.addWidget(self.label)
        layout.addStretch(1)
        layout.setSpacing(16)
        layout.setContentsMargins(0, 0, 0, 32)
        widget.setLayout(layout)

        return widget

    def updateLabels(self):
        path = self.machine().states.export.data.target
        count = self.machine().states.export.data.total
        self.confirm.setText(
            f'<b>Successfully exported {count} sequences to "{path.name}"</b>')

    def onEntry(self, event):
        super().onEntry(event)
        self.updateLabels()
