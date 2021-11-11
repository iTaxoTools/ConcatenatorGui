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

"""StepExport"""

from PySide6 import QtCore
from PySide6 import QtWidgets

from tempfile import TemporaryDirectory
from enum import Enum, IntEnum, auto
from datetime import datetime
from itertools import chain
from pathlib import Path

import shutil

from itaxotools.common.utility import AttrDict
from itaxotools.common.param.model import Model as ParamModel
from itaxotools.common.param.view import PlainView
from itaxotools.concatenator import (
    FileType, FileFormat, GeneStream, FileWriter,
    get_writer, get_extension, read_from_path)
from itaxotools.concatenator.library.operators import (
    OpChainGenes, OpTranslateGenes, OpApplyToGene, OpUpdateMetadata)

from .. import step_state_machine as ssm

from .wait import StepWaitBar


class FileScheme(Enum):
    InterNexus = ("Interleaved Nexus", FileType.File, FileFormat.Nexus)
    ConcatFasta = ("Concatenated Fasta", FileType.File, FileFormat.Fasta)
    ConcatPhylip = ("Concatenated Phylip", FileType.File, FileFormat.Phylip)
    ConcatAli = ("Concatenated Ali", FileType.File, FileFormat.Ali)
    MultiFasta = ("Multifile Fasta", FileType.Directory, FileFormat.Fasta)
    MultiPhylip = ("Multifile Phylip", FileType.Directory, FileFormat.Phylip)
    MultiAli = ("Multifile Ali", FileType.Directory, FileFormat.Ali)
    PartitionFinder = ("PartitionFinder", FileType.Directory, FileFormat.PartitionFinder)
    IQTree = ("IQTree", FileType.Directory, FileFormat.IQTree)
    ConcatTab = ("Tabfile", FileType.File, FileFormat.Tab)

    def __init__(self, text: str, type: FileType, format: FileFormat):
        self.text = text
        self.type = type
        self.format = format


class FileCompression(Enum):
    Uncompressed = ("None", FileType.File)
    ZipArchive = ("Zip Archive", FileType.ZipArchive)

    def __init__(self, text: str, type: FileType):
        self.text = text
        self.type = type


class CodonOptions(IntEnum):
    Mark = auto()
    Export = auto()


class OrderingOptions(IntEnum):
    SameAsInput = auto()
    Alphabetical = auto()


class StepExportEdit(ssm.StepTriStateEdit):

    description = 'Configure output'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.writer: Optional[FileWriter] = None
        self.scheme_changed()
        self.compression_changed()

    def onEntry(self, event):
        super().onEntry(event)
        self.footer.next.setText('&Export')

    def draw(self):
        widget = QtWidgets.QWidget()

        text = (
            'Please select the target file format '
            'and the desired options, then click "Export".')
        label = QtWidgets.QLabel(text)

        layout = QtWidgets.QGridLayout()
        layout.addWidget(label, 0, 0, 1, 4)
        layout.addLayout(self.draw_left(), 1, 0)
        layout.addLayout(self.draw_right(), 1, 1)
        layout.setColumnStretch(3, 1)
        layout.setRowStretch(3, 1)
        layout.setSpacing(16)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)

        return widget

    def draw_left(self):
        layout = QtWidgets.QGridLayout()

        scheme = QtWidgets.QComboBox()
        compression = QtWidgets.QComboBox()
        timestamp = QtWidgets.QCheckBox('Append timestamp to filename.')
        timestamp.setChecked(True)

        scheme.setToolTip((
            'Select one of the available sequence file formats.' + '\n'
            'Some options may only be available for certain formats.'
            ))
        compression.setToolTip((
            'If a compression type is selected, all output files' + '\n'
            'will be contained in a single archive of that type.'
            ))

        label_format = QtWidgets.QLabel('Output file format:')
        label_compression = QtWidgets.QLabel('File compression:')

        label_format.setToolTip(scheme.toolTip())
        label_compression.setToolTip(compression.toolTip())

        for item in FileScheme:
            scheme.addItem(item.text, item)
        for item in FileCompression:
            compression.addItem(item.text, item)

        scheme.currentIndexChanged.connect(self.scheme_changed)
        compression.currentIndexChanged.connect(self.compression_changed)

        layout.addWidget(label_format, 0, 0)
        layout.addWidget(label_compression, 1, 0)
        layout.addWidget(scheme, 0, 1)
        layout.addWidget(compression, 1, 1)
        layout.addWidget(timestamp, 3, 0, 1, 2)

        layout.setRowStretch(4, 1)
        layout.setColumnMinimumWidth(1, 160)
        layout.setSpacing(16)
        layout.setContentsMargins(24, 16, 24, 16)

        self.scheme = scheme
        self.compression = compression
        self.timestamp = timestamp

        return layout

    def draw_right(self):
        layout = QtWidgets.QVBoxLayout()

        self.param_view = PlainView()

        layout.addWidget(self.param_view)
        layout.addStretch(1)
        layout.setContentsMargins(24, 16, 24, 16)
        layout.setContentsMargins(24, 8, 24, 8)

        return layout

    def infer_writer(self):
        scheme = self.scheme.currentData()
        type = scheme.type
        compression = self.compression.currentData()
        if compression.type is not FileType.File:
            type = compression.type
        writer = get_writer(type, scheme.format)
        self.writer = writer
        if writer is None:
            raise Exception(f'Writer not found: {type} {scheme.format}')
        self.param_model = ParamModel(writer.params)
        self.param_view.setModel(self.param_model)

    def infer_dialog_filter(self):
        scheme = self.scheme.currentData()
        extension = get_extension(self.writer.type, self.writer.format)
        glob = f'*{extension}' if extension else '*'
        return f'{scheme.text} ({glob})'

    def infer_base_name(self):
        names = [path.stem for path in self.machine().states.input.data.files]
        if len(names) == 1:
            return names[0]
        return 'output'

    def get_time_string(self):
        if self.timestamp.isChecked():
            return '_' + datetime.now().strftime("%Y%m%dT%H%M%S")
        return ''

    def scheme_changed(self, index=0):
        scheme = self.scheme.currentData()
        scheme_is_dir = bool(scheme.type == FileType.Directory)
        self.compression.setEnabled(scheme_is_dir)
        if not scheme_is_dir:
            self.compression.setCurrentIndex(0)
        # Update option widgets here
        self.infer_writer()

    def compression_changed(self, index=0):
        # Update option widgets here
        self.infer_writer()


class StepExportWait(StepWaitBar):
    pass


class StepExportDone(ssm.StepTriStateDone):

    def onEntry(self, event):
        self.machine().navigate(ssm.NavigateAction.Next)


class StepExportFail(ssm.StepTriStateFail):
    description = 'Export failed'

    def onEntry(self, event):
        super().onEntry(event)
        message = f'{type(self.exception).__name__}'
        self.parent().update(text=f'Export failed: {message}')


class StepExport(ssm.StepTriState):
    title = 'Export sequence data'

    StepEdit = StepExportEdit
    StepWait = StepExportWait
    StepDone = StepExportDone
    StepFail = StepExportFail

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = AttrDict()
        self.data.temp = None
        self.data.target = None
        self.data.counter = 0
        self.data.total = 0

    def work(self):
        self.data.counter = 0
        self.data.total = len([
            cs for cs in self.machine().states.input.data.charsets.values()
            if cs.translation is not None])

        def checker_func(series):
            text = (
                f'Exporting sequence {self.data.counter}/{self.data.total}: '
                f'{series.name}')
            self.update(self.data.counter, self.data.total, text)
            with self.states.wait.redirect():
                print(text)
            QtCore.QThread.msleep(4)
            self.data.counter += 1
            self.worker.check()
            return series

        self.data.temp = TemporaryDirectory(prefix='concat_out_')
        out = Path(self.data.temp.name) / 'out'
        with self.states.wait.redirect():
            print('Exporting files, please wait...\n')
            print(f'Writing to temporary file: {out}\n')
        writer = self.states.edit.writer
        streams = [
            read_from_path(file.path)
            for file in self.machine().states.input.data.files.values()]
        cache = Path(self.machine().states.align_sets.temp_cache.name)
        all_streams = [read_from_path(cache)] + streams
        translation = self.machine().states.filter.translation
        stream = (
            GeneStream(chain(*all_streams))
            .pipe(OpChainGenes())
            .pipe(OpUpdateMetadata(self.machine().states.codons.metas))
            .pipe(OpTranslateGenes(translation))
            .pipe(OpApplyToGene(checker_func))
            )
        writer(stream, out)
        with self.states.wait.redirect():
            print(f'Done exporting, moving results to {self.data.target}')
        shutil.move(out, self.data.target)
        return self.data.total

    def onFail(self, exception, trace):
        self.states.wait.logio.writeline('')
        self.states.wait.logio.writeline(trace)
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Critical)
        msgBox.setText(type(exception).__name__)
        msgBox.setInformativeText(str(exception))
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.machine().parent().msgShow(msgBox)

    def filterNext(self, event):
        basename = self.states.edit.infer_base_name()
        basename += self.states.edit.get_time_string()
        (fileName, _) = QtWidgets.QFileDialog.getSaveFileName(
            self.machine().parent(),
            self.machine().parent().title + ' - Export',
            QtCore.QDir.currentPath() + '/' + basename,
            self.states.edit.infer_dialog_filter())
        if not fileName:
            return False
        if self.states.edit.writer.type == FileType.Directory:
            if Path(fileName).exists():
                msgBox = QtWidgets.QMessageBox(self.machine().parent())
                msgBox.setWindowTitle(self.machine().parent().title)
                msgBox.setIcon(QtWidgets.QMessageBox.Critical)
                msgBox.setText(
                    'Destination must be the name of a new directory')
                msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
                self.machine().parent().msgShow(msgBox)
                return False
        self.data.target = Path(fileName)
        return True

    def filterCancel(self, event):
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Question)
        msgBox.setText('Cancel file export?')
        msgBox.setStandardButtons(
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.No)
        res = self.machine().parent().msgShow(msgBox)
        return res == QtWidgets.QMessageBox.Yes
