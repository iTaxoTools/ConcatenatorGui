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

from tempfile import TemporaryDirectory
from itertools import chain
from typing import Optional
from pathlib import Path
from time import sleep

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from itaxotools.common.utility import AttrDict

from itaxotools.concatenator import (
    FileType, FileFormat, GeneStream, FileWriter,
    get_writer, get_extension, read_from_path)
from itaxotools.concatenator.library.model import GeneSeries, Operator
from itaxotools.concatenator.library.operators import OpUpdateMetadata
from itaxotools.concatenator.library.codons import (
    GeneticCode, ReadingFrame, _GC_DESCRIPTIONS)
from itaxotools.concatenator.library.operators import (
    OpChainGenes, OpTranslateGenes, OpApplyToGene, OpTagSet,
    OpUpdateMetadata, OpFilterGenes, OpGeneralInfo,
    OpGeneralInfoPerFile, OpGeneralInfoPerGene,
    OpGeneralInfoTagMafftRealigned,
    OpGeneralInfoTagPaddedLength,
    OpGeneralInfoTagPaddedCodonPosition)

from itaxotools.pygblocks import compute_mask, trim_sequence

from .. import model
from .. import widgets
from .. import step_state_machine as ssm
from ..trim import OpTrimGblocks, OpTrimClipKit

from .wait import StepWaitBar
from .align import RichRadioButton


class StepTrimEdit(ssm.StepTriStateEdit):

    description = 'Select trimming method'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = AttrDict()
        self.data.skip = False
        self.data.last = None

    def draw(self):
        widget = QtWidgets.QWidget()

        label = QtWidgets.QLabel('You may choose to trim your sequences. The same method will be used for all markers:')

        gblocks = RichRadioButton('Gblocks:', 'eliminate poorly aligned positions and divergent regions.', widget)
        clipkit = RichRadioButton('Clipkit:', 'only keep phylogenetically informative sites.', widget)
        skip = QtWidgets.QRadioButton('Skip trimming', widget)
        skip.setStyleSheet("QRadioButton { letter-spacing: 1px; }")
        clipkit.setChecked(True)

        radios = QtWidgets.QVBoxLayout()
        radios.addWidget(gblocks)
        radios.addWidget(clipkit)
        radios.addSpacing(16)
        radios.addWidget(skip)
        radios.setSpacing(16)
        radios.setContentsMargins(16, 0, 0, 0)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(label)
        layout.addLayout(radios)
        layout.addStretch(1)
        layout.setSpacing(24)
        layout.setContentsMargins(0, 0, 0, 32)
        widget.setLayout(layout)

        self.radios = AttrDict()
        self.radios.gblocks = gblocks
        self.radios.clipkit = clipkit
        self.radios.skip = skip

        return widget

    def get_method(self):
        for method in self.radios:
            if self.radios[method].isChecked():
                return method


class StepTrimWait(StepWaitBar):
    pass


class StepTrimDone(ssm.StepTriStateDone):
    description = 'Site trimming complete'

    def onEntry(self, event):
        super().onEntry(event)
        self.parent().update(
            text=f'Successfully trimmed sites.')


class StepTrimFail(ssm.StepTriStateFail):
    description = 'Position trimming failed'

    def onEntry(self, event):
        super().onEntry(event)
        message = f'{type(self.exception).__name__}'
        self.parent().update(
            text=f'Position trimming failed: {message}')


class StepTrim(ssm.StepTriState):
    title = 'Trim Sequence Sites'

    StepEdit = StepTrimEdit
    StepWait = StepTrimWait
    StepDone = StepTrimDone
    StepFail = StepTrimFail

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.temp_cache = None
        self.charsets_cached = set()
        self.method_last = None
        self.process = None

    def onEntry(self, event):
        super().onEntry(event)
        input_update = self.machine().states.input.timestamp_get()
        last_filter_update = self.machine().states.filter.timestamp_get()
        align_options_update = self.machine().states.align_options.timestamp_get()
        align_sets_update = self.machine().states.align_sets.timestamp_get()
        last_update = max(input_update, last_filter_update, align_options_update, align_sets_update)
        if last_update > self.timestamp_get():
            self.temp_cache = TemporaryDirectory(prefix='concat_trim_cache_')
            self.charsets_cached = set()
            self.method_last = None

    def onExit(self, event):
        super().onExit(event)
        method = self.states.edit.get_method()
        if method != self.method_last:
            self.timestamp_set()

    def work(self):
        charsets = {
            k for k, v in self.machine().states.input.data.charsets.items()
            if v.translation is not None and v.samples
            and v.name not in self.charsets_cached}
        method = self.states.edit.get_method()
        with self.states.wait.redirect():
            if not charsets and method == self.method_last:
                print('NOTHING TO DO')
                return
            stream = self.work_get_stream()
            self.work_trim(stream, method)
            self.method_last = method

    def work_trim(self, stream: GeneStream, method: str):
        self.update(0, 0, 'Getting ready...')
        print('Trimming sequences using pyGblocks:')
        print()
        print('Options:\n')
        print('- defaults')
        print(f'\n{"-"*20}\n')

        if method == 'gblocks':
            operator = OpTrimGblocks()
        elif method == 'clipkit':
            operator = OpTrimClipKit()
        else:
            raise Exception('Unexpected trimming method')

        stream = stream.pipe(operator)

        path = Path(self.temp_cache.name)
        writer = get_writer(FileType.Directory, FileFormat.Fasta)
        writer.params.translate_missing.value = ''
        writer.params.translate_gap.value = ''
        writer.params.padding.value = ''
        writer.params.sanitize.value = False
        writer.params.drop_empty.value = False
        writer(stream, path,)

        self.charsets_cached = {
            k for k, v in self.machine().states.input.data.charsets.items()
            if v.translation is not None and k in operator.genes}
        print('Done trimming!')
        print()
        self.update(1, 1, 'text')

        print(self.temp_cache.name)

    def work_get_stream(self):
        input_streams = [
            read_from_path(file.path)
            for file in self.machine().states.input.data.files.values()]
        aligned_cache = Path(self.machine().states.align_sets.temp_cache.name)
        aligned_charsets = {
            k for k, v in self.machine().states.input.data.charsets.items()
            if v.aligned and v.translation is not None}
        aligned_stream = (
            read_from_path(aligned_cache)
            .pipe(OpFilterGenes(aligned_charsets))
            )
        all_streams = [aligned_stream] + input_streams
        stream = (
            GeneStream(chain(*all_streams))
            .pipe(OpChainGenes())
            )
        return stream

    def skipWait(self):
        method = self.states.edit.get_method()
        return bool(method == 'skip')

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
