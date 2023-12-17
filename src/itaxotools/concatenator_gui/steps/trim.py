# -----------------------------------------------------------------------------
# ConcatenatorQt - GUI for Concatenator
# Copyright (C) 2021-2023  Patmanidis Stefanos
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

"""StepTrim"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui

from time import sleep

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from itaxotools.concatenator.library.model import GeneSeries
from itaxotools.concatenator.library.operators import OpUpdateMetadata
from itaxotools.concatenator.library.codons import (
    GeneticCode, ReadingFrame, _GC_DESCRIPTIONS)

from .. import model
from .. import widgets
from .. import step_state_machine as ssm

from .wait import StepWaitBar



class StepTrimEdit(ssm.StepTriStateEdit):

    description = 'Select trimming method'

    def onEntry(self, event):
        super().onEntry(event)
        last_filter_update = self.machine().states.filter.timestamp_get()
        if last_filter_update > self.timestamp_get():
            self.timestamp_set()

    def draw(self):
        widget = QtWidgets.QWidget()

        text = ('You may chose to trim your sequences using Gblocks.')
        label = QtWidgets.QLabel(text)

        frame = QtWidgets.QFrame()

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(label)
        layout.addWidget(frame, 1)
        layout.setSpacing(24)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)

        return widget


class StepTrimWait(StepWaitBar):
    pass


class StepTrimDone(ssm.StepTriStateDone):
    description = 'Position trimming complete'

    def onEntry(self, event):
        super().onEntry(event)
        self.parent().update(
            text=f'Successfully trimmed positions.')


class StepTrimFail(ssm.StepTriStateFail):
    description = 'Codon subsetting failed'

    def onEntry(self, event):
        super().onEntry(event)
        message = f'{type(self.exception).__name__}'
        self.parent().update(
            text=f'Position trimming failed: {message}')


class StepTrim(ssm.StepTriState):
    title = 'Trim Sequence Positions'

    StepEdit = StepTrimEdit
    StepWait = StepTrimWait
    StepDone = StepTrimDone
    StepFail = StepTrimFail

    def work(self):
        with self.states['wait'].redirect():
            self.update(0, 0, 'Getting ready...')
            print('Pretend to do some actual work...')
            sleep(2)
            print('Done pretending!')
            self.update(1, 1, 'text')

    def onFail(self, exception, trace):
        # raise exception
        self.states.wait.logio.writeline('')
        self.states.wait.logio.writeline(trace)
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Critical)
        msgBox.setText(type(exception).__name__)
        msgBox.setInformativeText(str(exception))
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.machine().parent().msgShow(msgBox)

    def filterCancel(self, event):
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Question)
        msgBox.setText('Cancel trimming?')
        msgBox.setStandardButtons(
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.No)
        res = self.machine().parent().msgShow(msgBox)
        return res == QtWidgets.QMessageBox.Yes
