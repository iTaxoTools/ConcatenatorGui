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

from itaxotools import common
from itaxotools.common.utility import AttrDict
from itaxotools.common.widgets import VLineSeparator, PushButton
from itaxotools.common.param.model import Model as ParamModel
from itaxotools.common.param.view import PlainView
from itaxotools.concatenator import (
    FileType, FileFormat, GeneStream, FileWriter,
    get_writer, get_extension, read_from_path)
from itaxotools.concatenator.library.operators import (
    OpChainGenes, OpTranslateGenes, OpApplyToGene, OpUpdateMetadata)
from itaxotools import mafftpy
from itaxotools import fasttreepy
from itaxotools.fasttreepy.gui.main import CustomView as TreeParamView
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


class TreeOptionsDialog(QtWidgets.QDialog):
    def __init__(self, params, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.params = params
        self.draw()

        self.setWindowTitle(self.parent().title)
        # self.setFixedSize(self.sizeHint())

    def draw(self):
        self.model = ParamModel(self.params)
        self.view = TreeParamView(self.model, showResetButton=False)
        self.view.widget().layout().setContentsMargins(16, 12, 8, 12)
        self.view.setStyleSheet("""
            QScrollArea {
                border: 0px;
                border-bottom: 1px solid Palette(Mid);
                }
        """)

        ok = common.widgets.PushButton('OK')
        ok.clicked.connect(self.accept)
        ok.setDefault(True)
        reset = common.widgets.PushButton('Reset')
        reset.clicked.connect(self.reset)

        buttons = QtWidgets.QHBoxLayout()
        buttons.addSpacing(16)
        buttons.addWidget(reset)
        buttons.addWidget(ok)
        buttons.addSpacing(16)
        buttons.setSpacing(8)
        buttons.setContentsMargins(0, 0, 0, 0)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.view)
        layout.addSpacing(8)
        layout.addLayout(buttons)
        layout.addSpacing(16)
        layout.setContentsMargins(0, 0, 0, 0)
        # layout.setContentsMargins(24, 16, 24, 16)
        self.setLayout(layout)

    def reset(self):
        self.model.resetParams()


class StepExportEdit(ssm.StepTriStateEdit):

    description = 'Configure output'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.phylo_params = fasttreepy.params.params()
        self.phylo_available = False
        self.writer: Optional[FileWriter] = None
        self.scheme_changed()
        self.compression_changed()

    def onEntry(self, event):
        super().onEntry(event)
        self.footer.next.setText('&Export')
        if self.checkGenesChanged():
            self.updatePhyloAvailable()
            self.updatePhyloLayout()
            self.handlePhyloUpdate()

    def draw(self):
        widget = QtWidgets.QWidget()

        text = (
            'Export using the desired file format. '
            'Hover options for more information.')
        label = QtWidgets.QLabel(text)
        text_trees = (
            'You may additionally calculate phylogenetic trees '
            'using FastTree.')
        label_trees = QtWidgets.QLabel(text_trees)
        separator = VLineSeparator(1)

        layout = QtWidgets.QGridLayout()
        layout.addWidget(label, 0, 0)
        layout.addWidget(separator, 0, 2, 6, 1)
        layout.addLayout(self.draw_left(), 1, 0)
        layout.addLayout(self.draw_right(), 1, 3, 5, 1)
        layout.addWidget(label_trees, 3, 0)
        layout.addLayout(self.draw_trees(), 4, 0)
        layout.setColumnStretch(1, 1)
        layout.setRowMinimumHeight(2, 8)
        layout.setRowStretch(5, 1)
        layout.setSpacing(8)
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
        timestamp.setToolTip(
            'The timestamp follows ISO 8601 in local time.')

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
        self.param_view.layout().setSpacing(16)
        self.param_view.layout().setContentsMargins(8, 16, 0, 0)

        layout.addWidget(self.param_view)
        layout.addStretch(1)
        layout.setContentsMargins(0, 0, 0, 0)

        return layout

    def draw_trees(self):
        layout = QtWidgets.QVBoxLayout()

        phylo_warning = QtWidgets.QLabel(
            '\u21AA Only available if all sequences are aligned/uniform.')

        phylo_concat = QtWidgets.QCheckBox(
            'Calculate tree for the concatenated alignment.')

        phylo_all = QtWidgets.QCheckBox(
            'Calculate trees for each single-gene alignment.')

        text = 'Requires that all genes are aligned (eg. with MAFFT).'
        phylo_concat.setToolTip(text)
        phylo_all.setToolTip(text)

        phylo_concat.stateChanged.connect(self.handlePhyloUpdate)
        phylo_all.stateChanged.connect(self.handlePhyloUpdate)

        phylo_config = PushButton(
            'FastTree Options', onclick=self.handlePhyloShowConfigDialog)
        phylo_config.setMaximumWidth(160)

        layout.addWidget(phylo_warning)
        layout.addWidget(phylo_concat)
        layout.addWidget(phylo_all)
        layout.addSpacing(8)
        layout.addWidget(phylo_config)
        layout.setSpacing(16)
        layout.setContentsMargins(24, 12, 24, 16)

        self.phylo_warning = phylo_warning
        self.phylo_concat = phylo_concat
        self.phylo_all = phylo_all
        self.phylo_config = phylo_config

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
        self.infer_writer()

    def compression_changed(self, index=0):
        self.infer_writer()

    def checkGenesChanged(self):
        last_update = max(
            self.machine().states.filter.timestamp_get(),
            self.machine().states.align_sets.timestamp_get())
        return self.timestamp_get() < last_update

    def updatePhyloAvailable(self):
        self.phylo_available = not any(
            cs for cs in self.machine().states.input.data.charsets.values()
            if cs.translation is not None and not (
                cs.aligned or cs.uniform == 'Yes'
            ))
        self.timestamp_set()

    def updatePhyloLayout(self):
        enabled = self.phylo_available
        self.phylo_config.setVisible(enabled)
        self.phylo_warning.setVisible(not enabled)
        self.phylo_concat.setVisible(enabled)
        self.phylo_all.setVisible(enabled)
        if not enabled:
            self.phylo_concat.setChecked(False)
            self.phylo_all.setChecked(False)

    def handlePhyloUpdate(self):
        self.phylo_config.setEnabled(
            self.phylo_concat.isChecked() or
            self.phylo_all.isChecked())

    def handlePhyloShowConfigDialog(self):
        self.dialog = TreeOptionsDialog(
            self.phylo_params, self.machine().parent())
        self.dialog.setModal(True)
        self.dialog.show()


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
        self.data.seqs = 0
        self.data.trees = 0
        self.data.phylo_prep = None

    def work(self):
        self.data.total = 0
        self.data.trees = 0
        self.data.phylo_do_concat = self.states.edit.phylo_concat.isChecked()
        self.data.phylo_do_all = self.states.edit.phylo_all.isChecked()
        # self.data.phylo_do_concat = True  # dev
        # self.data.phylo_do_all = True  # dev
        if self.data.phylo_do_concat or self.data.phylo_do_all:
            self.work_export('Step 1/3: Exporting sequences')
            self.work_phylo_prep('Step 2/3: Phylogeny preparation')
            self.work_phylo_calc('Step 3/3: Phylogeny calculation')
        else:
            self.work_export('Please wait...')

    def checker_func(self, series):
        text = (
            f'Exporting sequence {self.data.counter}/{self.data.total}: '
            f'{series.name}')
        self.update(self.data.counter, self.data.total, text)
        with self.states.wait.redirect():
            print(text)
        # QtCore.QThread.msleep(400)
        self.data.counter += 1
        self.worker.check()
        return series

    def work_get_stream(self):
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
            .pipe(OpApplyToGene(self.checker_func))
            )
        return stream

    def work_export(self, description):
        self.header.showTask(title=self.title, description=description)
        self.data.counter = 0
        self.data.total = len([
            cs for cs in self.machine().states.input.data.charsets.values()
            if cs.translation is not None])
        self.data.temp = TemporaryDirectory(prefix='concat_out_')
        out = Path(self.data.temp.name) / 'out'
        with self.states.wait.redirect():
            print('Exporting files, please wait...\n')
            print(f'Writing to temporary file: {out}\n')
        stream = self.work_get_stream()
        writer = self.states.edit.writer
        writer(stream, out)
        # QtCore.QThread.msleep(200)
        with self.states.wait.redirect():
            print(f'Done exporting, moving results to {self.data.target}')
        shutil.move(out, self.data.target)
        self.data.seqs = self.data.total

    def work_phylo_prep(self, description):
        self.header.showTask(title=self.title, description=description)
        if self.data.phylo_do_concat:
            self.data.trees += 1
        if self.data.phylo_do_all:
            self.data.trees += self.data.total
        self.data.total = self.data.seqs * sum([
            self.data.phylo_do_concat, self.data.phylo_do_all])
        self.data.counter = 0
        self.data.phylo_prep = TemporaryDirectory(prefix='concat_phylo_prep_')
        self.work_phylo_prep_concat()
        self.work_phylo_prep_all()

    def work_phylo_prep_concat(self):
        out = Path(self.data.phylo_prep.name) / 'concat'
        with self.states.wait.redirect():
            print('\nPreparing for concatenated phylogeny, please wait...\n')
            print(f'Writing to temporary file: {out}\n')
        stream = self.work_get_stream()
        writer = get_writer(FileType.File, FileFormat.Fasta)
        writer(stream, out)
        # QtCore.QThread.msleep(200)
        return 1

    def work_phylo_prep_all(self):
        out = Path(self.data.phylo_prep.name) / 'all'
        with self.states.wait.redirect():
            print('\nPreparing for individual phylogeny, please wait...\n')
            print(f'Writing to temporary file: {out}\n')
        stream = self.work_get_stream()
        writer = get_writer(FileType.Directory, FileFormat.Fasta)
        writer(stream, out)
        # QtCore.QThread.msleep(200)
        return self.data.total

    def work_phylo_calc(self, description):
        self.header.showTask(title=self.title, description=description)
        # QtCore.QThread.msleep(200)

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
