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
from PySide6 import QtGui

from tempfile import TemporaryDirectory
from enum import Enum, IntEnum, auto
from datetime import datetime
from itertools import chain
from pathlib import Path
from dataclasses import dataclass

import shutil
import pandas as pd

from itaxotools import common
from itaxotools.common.utility import AttrDict
from itaxotools.common.widgets import VLineSeparator, PushButton
from itaxotools.common.param.model import Model as ParamModel
from itaxotools.common.param.view import PlainView
from itaxotools.concatenator import (
    FileType, FileFormat, GeneStream, FileWriter,
    get_writer, get_extension, read_from_path)
from itaxotools.concatenator.library.file_utils import ZipFile, ZipPath
from itaxotools.concatenator.library.operators import (
    OpChainGenes, OpTranslateGenes, OpApplyToGene, OpTagSet,
    OpUpdateMetadata, OpFilterGenes, OpGeneralInfo,
    OpGeneralInfoPerFile, OpGeneralInfoPerGene,
    OpGeneralInfoTagMafftRealigned,
    OpGeneralInfoTagPaddedLength,
    OpGeneralInfoTagPaddedCodonPosition)
from itaxotools import mafftpy
from itaxotools.fasttreepy import PhylogenyApproximation
from itaxotools.fasttreepy.params import params as fasttreepy_params
from itaxotools.fasttreepy.gui.main import CustomView as TreeParamView
from .. import step_state_machine as ssm

from .wait import StepWaitBar


def work_fasttree(src, out, param):
    from sys import stdout
    phylo = PhylogenyApproximation(str(src))
    phylo.param = param
    phylo.target = out
    phylo.log = stdout
    phylo.fetch = lambda: out  # ugly! fix api instead
    phylo.run()


def copy_file_to_archive(src, dst):
    with src.open('r') as fin, dst.open('w') as fout:
        for line in fin:
            fout.write(line)


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


@dataclass
class DiagnoserParams:
    report: bool = True
    disjoint: bool = True
    foreign: bool = True
    outliers: bool = False
    iqr: float = 20.0


class Diagnoser:
    def __init__(self, params: DiagnoserParams):
        self.params = params
        self.op_general_info = OpGeneralInfo()
        self.op_general_info_per_gene = OpGeneralInfoPerGene()
        self.op_general_info_per_file = OpGeneralInfoPerFile()
        self.temp = TemporaryDirectory(prefix='concat_diagnose_')
        self.filename = None
        self.timestamp = None
        self.scheme = None
        self.input_files_count = None
        self.disjoint_groups = None
        self.foreign_pairs = None
        self._writer_padding = False
        self._writer_frames = False

    def pipe_input_streams(self, input_streams):
        self.input_files_count = len(input_streams)
        if not self.params.report:
            return input_streams
        return [
            stream.pipe(self.op_general_info_per_file)
            for stream in input_streams]

    def pipe_aligned_stream(self, aligned_stream):
        if not self.params.report:
            return aligned_stream
        return aligned_stream.pipe(OpGeneralInfoTagMafftRealigned())

    def _pipe_general_info(self, stream):
        if not self.params.report and not self.params.disjoint:
            return stream
        return (
            stream
            .pipe(self.op_general_info)
            .pipe(self.op_general_info_per_gene)
            )

    def _pipe_padding_info(self, stream):
        if not self.params.report:
            return stream
        if self._writer_padding:
            stream = stream.pipe(OpGeneralInfoTagPaddedLength())
        if self._writer_frames:
            stream = stream.pipe(OpGeneralInfoTagPaddedCodonPosition())
        return stream

    def modify_writer_filters(self, writer):
        self._writer_padding = bool(getattr(writer, 'padding', False))
        self._writer_frames = bool(getattr(writer, 'adjust_frames', False))
        writer.filters.insert(0, self._pipe_padding_info)
        writer.filters.append(self._pipe_general_info)

    def get_table_header(self):
        return pd.Series({
            'Name of output file': self.filename,
            'Timestamp of output file': self.timestamp,
            'Output file format': self.scheme.text,
            'Number of input files used': self.input_files_count,
        })

    def get_dataframe(self):
        return self.op_general_info.table.dataframe

    def get_table_total(self):
        header = self.get_table_header()
        table = self.op_general_info.table.total_data()
        table["GC content of sequences"] = f'{table["GC content of sequences"]:.1f}'
        table["% missing data"] = f'{table["% missing data"]:.1f}'
        return pd.concat([header, table])

    def get_table_by_taxon(self):
        return self.op_general_info.table.by_taxon()

    def get_table_by_gene(self):
        table = self.op_general_info.table
        return self.op_general_info_per_gene.get_info(table)

    def get_table_by_input_file(self):
        return self.op_general_info_per_file.get_info()

    def export_dir(self):
        return Path(self.temp.name)

    def export_table(
        self, table, name, header=True, index=True,
        sep='\t', float_format='%.1f'
    ):
        table.to_csv(
            self.export_dir() / name, header=header, index=index,
            sep=sep, float_format=float_format)

    def export_disjoint(self, name):
        self.disjoint_groups = 0
        with open(self.export_dir() / name, 'w') as file:
            for group in self.op_general_info.table.disjoint_taxon_groups():
                self.disjoint_groups += 1
                print(f'Group {self.disjoint_groups}', file=file)
                print('-------------', file=file)
                for taxon in group:
                    print(taxon, file=file)
                print('', file=file)

    def export_foreign(self, name):
        self.foreign_pairs = 0
        with open(self.export_dir() / name, 'w') as file:
            for x, y in self.op_general_info.table.unconnected_taxons():
                self.foreign_pairs += 1
                print(f'{x}\t{y}', file=file)

    def export_all(self):
        temp = self.export_dir()
        if self.params.report:
            self.export_table(self.get_table_total(), 'total_data.tsv', header=False)
            self.export_table(self.get_table_by_input_file(), 'by_input_file.tsv', index=False)
            self.export_table(self.get_table_by_taxon(), 'by_taxon.tsv')
            self.export_table(self.get_table_by_gene(), 'by_gene.tsv')
        if self.params.disjoint:
            self.export_disjoint('disjoint_groups.txt')
        if self.params.foreign:
            self.export_foreign('foreign_pairs.txt')


class DiagnoserOptionsDialog(QtWidgets.QDialog):
    def __init__(self, params: DiagnoserParams, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.params = params
        self.draw()
        self.load_params()

        self.setWindowTitle(self.parent().title)
        self.setFixedSize(self.sizeHint())

    def draw(self):
        self.header = QtWidgets.QLabel(
            'Diagnostic options for data validation: ')

        self.report = QtWidgets.QCheckBox(
            'Produce summary report tables.')
        self.disjoint = QtWidgets.QCheckBox(
            'Detect disjoint sample groups.')
        self.foreign = QtWidgets.QCheckBox(
            'Detect foreign samples.')

        self.outliers = QtWidgets.QCheckBox(
            'Detect outlier sequences.')
        self.outliers.stateChanged.connect(self.setIqrEnabled)

        self.iqr_label = QtWidgets.QLabel(
            ' \u21AA IQR coefficient:   ')
        self.iqr_edit = QtWidgets.QLineEdit('1.0')
        validator = QtGui.QDoubleValidator(self.iqr_edit)
        validator.setBottom(0)
        self.iqr_edit.setValidator(validator)

        self.footer = QtWidgets.QLabel(
            'Hover fields for more details.')

        self.report.setToolTip(
            'Produce summary tables per sample, marker and input file. \n'
            'Also produces a short summary for the whole dataset.')
        self.disjoint.setToolTip(
            'List groups of samples, such that samples in different groups share no markers. \n'
            'A warning will be issued if there are more than one such groups.')
        self.foreign.setToolTip(
            'Find pairs of samples that have no markers in common. \n'
            'A warning will be issued if any exist.')
        self.outliers.setToolTip(
            'List outlier samples using SequenceBouncer. \n'
            'A warning will be issued if any exist.')

        iqr_tip = (
            'Coefficient multiplied by the interquartile range that helps \n'
            'define an outlier sequence (default is 1.0).')
        self.iqr_label.setToolTip(iqr_tip)
        self.iqr_edit.setToolTip(iqr_tip)

        ok = common.widgets.PushButton('OK')
        ok.clicked.connect(self.accept)
        ok.setDefault(True)
        cancel = common.widgets.PushButton('Cancel')
        cancel.clicked.connect(self.reject)

        buttons = QtWidgets.QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(cancel)
        buttons.addWidget(ok)
        buttons.setSpacing(8)
        buttons.setContentsMargins(0, 0, 0, 0)

        options = QtWidgets.QGridLayout()
        options.setRowMinimumHeight(10, 8)
        options.addWidget(self.header, 11, 0, 1, 4)
        options.setRowMinimumHeight(20, 16)
        options.addWidget(self.report, 21, 1, 1, 3)
        options.addWidget(self.disjoint, 22, 1, 1, 3)
        options.addWidget(self.foreign, 23, 1, 1, 3)
        options.setRowMinimumHeight(30, 16)
        options.addWidget(self.outliers, 31, 1, 1, 3)
        options.addWidget(self.iqr_label, 32, 1, 1, 1)
        options.addWidget(self.iqr_edit, 32, 2, 1, 1)
        options.setRowMinimumHeight(40, 16)
        options.addWidget(self.footer, 41, 0, 1, 4)
        options.setRowMinimumHeight(50, 24)

        options.setColumnStretch(3, 1)
        options.setColumnMinimumWidth(0, 16)
        options.setColumnMinimumWidth(4, 16)
        options.setContentsMargins(0, 0, 0, 0)

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(options)
        layout.addLayout(buttons)
        layout.setContentsMargins(24, 16, 24, 16)
        self.setLayout(layout)

    def setIqrEnabled(self, state):
        self.iqr_label.setEnabled(state)
        self.iqr_edit.setEnabled(state)

    def accept(self):
        super().accept()
        self.save_params()

    def load_params(self):
        self.report.setChecked(self.params.report)
        self.disjoint.setChecked(self.params.disjoint)
        self.foreign.setChecked(self.params.foreign)
        self.outliers.setChecked(self.params.outliers)
        self.iqr_edit.setText(str(self.params.iqr))
        self.setIqrEnabled(self.params.outliers)

    def save_params(self):
        self.params.report = self.report.isChecked()
        self.params.disjoint = self.disjoint.isChecked()
        self.params.foreign = self.foreign.isChecked()
        self.params.outliers = self.outliers.isChecked()
        self.params.iqr = float(self.iqr_edit.text())


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
        self.phylo_params = fasttreepy_params()
        self.diagnoser_params = DiagnoserParams()
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
        diagnose_config = PushButton(
            'Diagnostic Options', onclick=self.handleDiagnoseShowConfigDialog)
        diagnose_config.setMaximumWidth(160)

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
        diagnose_config.setToolTip(
            'Select between a series of diagnostic options'
            'for data validation.')

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
        layout.addWidget(diagnose_config, 5, 0, 1, 2)

        layout.setRowStretch(5, 1)
        layout.setRowMinimumHeight(4, 8)
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
        return 'sequences'

    def get_timestamp(self):
        return datetime.now().strftime("%Y%m%dT%H%M%S")

    def get_time_string(self, timestamp=None):
        timestamp = timestamp or self.get_timestamp()
        if self.timestamp.isChecked():
            return '_' + timestamp
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

    def handleDiagnoseShowConfigDialog(self):
        self.dialog = DiagnoserOptionsDialog(
            self.diagnoser_params, self.machine().parent())
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
    title = 'Export Sequence Data'

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
        self.data.phylo_calc = None

    def work(self):
        self.data.total = 0
        self.data.trees = 0
        if self.data.phylo_do_concat or self.data.phylo_do_all:
            self.work_export('Step 1/3: Exporting sequences')
            self.work_phylo_prep('Step 2/3: Phylogeny preparation')
            self.work_phylo_calc('Step 3/3: Phylogeny calculation')
        else:
            self.work_export('Please wait...')
        self.work_save()

    def checker_func(self, series):
        text = (
            f'Exporting sequence {self.data.counter}/{self.data.total}: '
            f'{series.name}')
        self.update(self.data.counter, self.data.total, text)
        with self.states.wait.redirect():
            print(text)
        self.data.counter += 1
        self.worker.check()
        return series

    def work_get_stream(self):
        input_streams = [
            read_from_path(file.path)
            for file in self.machine().states.input.data.files.values()]
        input_streams = self.data.diagnoser.pipe_input_streams(input_streams)
        aligned_cache = Path(self.machine().states.align_sets.temp_cache.name)
        aligned_charsets = {
            k for k, v in self.machine().states.input.data.charsets.items()
            if v.aligned and v.translation is not None}
        aligned_stream = (
            read_from_path(aligned_cache)
            .pipe(OpFilterGenes(aligned_charsets))
            )
        aligned_stream = self.data.diagnoser.pipe_aligned_stream(aligned_stream)
        all_streams = [aligned_stream] + input_streams
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
        self.data.diagnoser.modify_writer_filters(writer)
        writer(stream, out)
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
        if self.data.phylo_do_concat:
            self.work_phylo_prep_concat()
        if self.data.phylo_do_all:
            self.work_phylo_prep_all()

    def work_phylo_prep_concat(self):
        out = Path(self.data.phylo_prep.name) / 'concat'
        with self.states.wait.redirect():
            print('\nPreparing for concatenated phylogeny, please wait...\n')
            print(f'Writing to temporary file: {out}\n')
        stream = self.work_get_stream()
        writer = get_writer(FileType.File, FileFormat.Fasta)
        writer(stream, out)
        return 1

    def work_phylo_prep_all(self):
        out = Path(self.data.phylo_prep.name) / 'all'
        with self.states.wait.redirect():
            print('\nPreparing for individual phylogeny, please wait...\n')
            print(f'Writing to temporary file: {out}\n')
        stream = self.work_get_stream()
        writer = get_writer(FileType.Directory, FileFormat.Fasta)
        writer(stream, out)
        return self.data.total

    def work_phylo_fasttree(self, src, out):
        # print('FASTTREE', src, '->', out)
        loop = QtCore.QEventLoop()
        self.process = common.threading.Process(
            work_fasttree, src, out, self.states.edit.phylo_params)
        self.process.setStream(self.states.wait.logio)
        self.process.done.connect(loop.quit)
        # process.fail does not send out the trace
        # self.process.fail.connect(self.worker.fail)
        self.process.fail.connect(
            lambda e: self.worker.fail.emit(e, '???'))
        self.process.start()
        loop.exec()
        self.worker.check()

    def work_phylo_calc(self, description):
        self.header.showTask(title=self.title, description=description)
        self.data.phylo_calc = TemporaryDirectory(prefix='concat_phylo_calc_')
        if self.data.phylo_do_concat:
            src = Path(self.data.phylo_prep.name) / 'concat'
            out = Path(self.data.phylo_calc.name) / 'concat'
            with self.states.wait.redirect():
                print('\nCalculating concatenated phylogeny...\n')
            self.work_phylo_fasttree(src, out)
        if self.data.phylo_do_all:
            src = Path(self.data.phylo_prep.name) / 'all'
            out = Path(self.data.phylo_calc.name) / 'all'
            out.mkdir()
            for item in src.iterdir():
                self.states.wait.reset()
                with self.states.wait.redirect():
                    print(f'\nCalculating phylogeny for {item.stem}...\n')
                self.work_phylo_fasttree(
                    src / item.name,
                    out / f'{item.stem}.tre')

    def work_save(self):
        self.states.wait.reset()
        with self.states.wait.redirect():
            print(f'Saving diagnostics...')

        self.data.diagnoser.export_all()

        with self.states.wait.redirect():
            print(f'Done exporting, moving results to {self.data.target}')

        target_out = self.data.target
        do_phylo = self.data.phylo_do_concat or self.data.phylo_do_all
        type = self.states.edit.writer.type
        if do_phylo and type == FileType.File:
            target_out.mkdir()
            basename = self.states.edit.infer_base_name()
            ext = self.states.edit.writer.format.extension
            target_out = target_out / f'{basename}{ext}'
        out = Path(self.data.temp.name) / 'out'
        shutil.move(out, target_out)

        if do_phylo:
            if type == FileType.File:
                target_phylo = target_out.parent
                move_func = shutil.move
            elif type == FileType.Directory:
                target_phylo = target_out
                move_func = shutil.move
            elif type == FileType.ZipArchive:
                archive = ZipFile(target_out, 'a')
                target_phylo = ZipPath(archive)
                move_func = copy_file_to_archive
        if self.data.phylo_do_concat:
            out = Path(self.data.phylo_calc.name) / 'concat'
            basename = self.states.edit.infer_base_name()
            move_func(out, target_phylo / f'{basename}.tre')
        if self.data.phylo_do_all:
            out = Path(self.data.phylo_calc.name) / 'all'
            for item in out.iterdir():
                move_func(item, target_phylo / item.name)

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

    def isTargetDir(self):
        type = self.states.edit.writer.type
        return (
            type == FileType.Directory or
            type == FileType.File and (
                self.data.phylo_do_concat or self.data.phylo_do_all))

    def filterNext(self, event):
        self.data.phylo_do_concat = self.states.edit.phylo_concat.isChecked()
        self.data.phylo_do_all = self.states.edit.phylo_all.isChecked()
        timestamp = self.states.edit.get_timestamp()
        basename = self.states.edit.infer_base_name()
        basename += self.states.edit.get_time_string(timestamp)
        (fileName, _) = QtWidgets.QFileDialog.getSaveFileName(
            self.machine().parent(),
            self.machine().parent().title + ' - Export',
            QtCore.QDir.currentPath() + '/' + basename,
            self.states.edit.infer_dialog_filter())
        if not fileName:
            return False
        if self.isTargetDir():
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
        self.data.diagnoser = Diagnoser(self.states.edit.diagnoser_params)
        self.data.diagnoser.scheme = self.states.edit.scheme.currentData()
        self.data.diagnoser.timestamp = timestamp
        self.data.diagnoser.filename = self.data.target.name
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
