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

"""StepWaitBar"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui

from contextlib import contextmanager

import sys

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from .. import step_state_machine as ssm


class ProgressBar(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.label = QtWidgets.QLabel('No task')
        self.label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        self.bar = QtWidgets.QProgressBar()
        self.bar.setTextVisible(False)
        self.bar.setStyleSheet("""
            QProgressBar {
                background: Palette(Light);
                border: 1px solid Palette(Mid);
                border-radius: 0px;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: Palette(Highlight);
                width: 1px;
            }""")
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.bar)
        layout.setSpacing(6)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)


class StepWaitBar(ssm.StepTriStateWait):
    description = 'Please wait...'

    def draw(self):
        widget = QtWidgets.QWidget()

        self.progress = ProgressBar()

        self.logger = common.widgets.TextEditLogger()
        self.logger.setFont(QtGui.QFontDatabase.systemFont(
            QtGui.QFontDatabase.FixedFont))
        self.logger.document().setDocumentMargin(10)
        self.logio = common.io.TextEditLoggerIO(self.logger)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.progress)
        layout.addWidget(self.logger, 1)
        layout.setSpacing(16)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)

        return widget

    def onEntry(self, event):
        super().onEntry(event)
        self.logger.clear()

    def update(self, val=None, max=None, text=None):
        if max is not None:
            self.progress.bar.setMaximum(max)
        if val is not None:
            self.progress.bar.setValue(val)
        if text is not None:
            self.progress.label.setText(text)

    @contextmanager
    def redirect(self):
        with common.io.redirect(sys, 'stdout', self.logio) as out:
            yield out
