# -----------------------------------------------------------------------------
# ConcatenatorQt - GUI for Concatenator
# Copyright (C) 2021-2024  Patmanidis Stefanos
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

"""StepAbout"""

from PySide6 import QtWidgets

from itaxotools import common
import itaxotools.common.resources # noqa

from .. import step_state_machine as ssm
from .. import widgets


class StepAbout(ssm.StepState):

    def draw(self):
        widget = QtWidgets.QWidget()

        path = common.resources.get(
            'itaxotools.concatenator_gui', 'docs/about.html')
        self.about = widgets.HtmlLabel(path)

        path = common.resources.get(
            'itaxotools.concatenator_gui', 'docs/cite.html')
        self.cite = widgets.HtmlLabel(path)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.about)
        layout.addStretch(1)
        layout.addWidget(self.cite)
        layout.setSpacing(24)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)

        return widget
