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

from enum import Enum, IntEnum, auto
from pathlib import Path

from itaxotools.common.utility import AttrDict
from itaxotools.concatenator import (
    FileType, FileFormat, get_writer, get_extension)

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

    def onEntry(self, event):
        super().onEntry(event)
        self.footer.next.setText('&Export')
        self.format_changed()

    def draw(self):
        widget = QtWidgets.QWidget()

        text = (
            'Please select the output file format and compression type, '
            'then click "Export" to save results:')
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
        order = QtWidgets.QComboBox()
        compression = QtWidgets.QComboBox()

        scheme.setToolTip((
            'Select one of the available sequence file formats.' + '\n'
            'Some options may only be available for certain formats.'
            ))
        order.setToolTip((
            'The order in which the character sets are written' + '\n'
            'in the output file information blocks.'
            ))
        compression.setToolTip((
            'If a compression type is selected, all output files' + '\n'
            'will be contained in a single archive of that type.'
            ))

        label_format = QtWidgets.QLabel('Output file format:')
        label_compression = QtWidgets.QLabel('File compression:')
        label_order = QtWidgets.QLabel('Character set order:')

        label_format.setToolTip(scheme.toolTip())
        label_compression.setToolTip(compression.toolTip())
        label_order.setToolTip(order.toolTip())

        for item in FileScheme:
            scheme.addItem(item.text, item)
        for item in FileCompression:
            compression.addItem(item.text, item)
        order.addItem('Alphabetical', OrderingOptions.Alphabetical)
        order.addItem('Same as input', OrderingOptions.SameAsInput)

        scheme.activated.connect(self.format_changed)
        compression.activated.connect(self.format_changed)

        label_order.setVisible(False)
        order.setVisible(False)

        layout.addWidget(label_format, 0, 0)
        layout.addWidget(label_compression, 1, 0)
        layout.addWidget(label_order, 2, 0)
        layout.addWidget(scheme, 0, 1)
        layout.addWidget(compression, 1, 1)
        layout.addWidget(order, 2, 1)

        layout.setRowStretch(3, 1)
        layout.setColumnMinimumWidth(1, 160)
        layout.setSpacing(16)
        layout.setContentsMargins(24, 16, 24, 16)

        self.scheme = scheme
        self.order = order
        self.compression = compression

        return layout

    def draw_right(self):
        layout = QtWidgets.QVBoxLayout()

        mark = QtWidgets.QRadioButton(
            'Write codon sets within output information block.')
        mark.setToolTip((
            'Codon sets will be written within:' + '\n'
            '- the sets block for Interleaved Nexus' + '\n'
            '- the cfg file for Partitionfinder'
            ))

        export = QtWidgets.QRadioButton(
            'Export codon sets as new character sets.')
        export.setToolTip(
            'Codon sets will be treated as new character sets.'
            )

        options = QtWidgets.QButtonGroup()
        options.addButton(mark)
        options.setId(mark, CodonOptions.Mark)
        options.addButton(export)
        options.setId(export, CodonOptions.Export)

        exclude = QtWidgets.QCheckBox(
            'Exclude sequences with missing data.')
        exclude.setToolTip((
            'Samples that contain only { N ? - } will not' + '\n'
            'be included in the respective output files.'
            ))

        mark.setChecked(True)
        exclude.setChecked(True)

        mark.setVisible(False)
        export.setVisible(False)
        exclude.setVisible(False)

        layout.addWidget(mark)
        layout.addWidget(export)
        layout.addSpacing(16)
        layout.addWidget(exclude)
        layout.addStretch(1)
        layout.setSpacing(8)
        layout.setContentsMargins(24, 16, 24, 16)

        self.codonOptions = options
        self.exclude = exclude

        return layout

    def infer_writer(self):
        scheme = self.scheme.currentData()
        type = scheme.type
        compression = self.compression.currentData()
        if compression.type is not FileType.File:
            type = compression.type
        return get_writer(type, scheme.format)

    def infer_dialog_filter(self):
        writer = self.infer_writer()
        scheme = self.scheme.currentData()
        extension = get_extension(writer.type, writer.format)
        glob = f'*{extension}' if extension else '*'
        return f'{scheme.text} ({glob})'

    def format_changed(self, index=0):
        scheme = self.scheme.currentData()
        allow_compression = bool(scheme.type == FileType.Directory)
        self.compression.setEnabled(allow_compression)
        if not allow_compression:
            self.compression.setCurrentIndex(0)
        # Update option widgets here


class StepExportWait(StepWaitBar):
    pass


class StepExportDone(ssm.StepTriStateDone):

    def onEntry(self, event):
        self.machine().navigate(ssm.NavigateAction.Next)


class StepExportFail(ssm.StepTriStateFail):
    description = 'Export failed'

    def onEntry(self, event):
        super().onEntry(event)
        self.parent().update(
            text=f'Export failed: {str(self.exception)}')


class StepExport(ssm.StepTriState):
    title = 'Export sequence data'

    StepEdit = StepExportEdit
    StepWait = StepExportWait
    StepDone = StepExportDone
    StepFail = StepExportFail

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = AttrDict()
        self.data.target = None
        self.data.count = 0

    def work(self):
        with self.states.wait.redirect():
            print('')
            return 0

    def filterNext(self, event):
        if self.states.edit.infer_writer().type == FileType.Directory:
            fileName = QtWidgets.QFileDialog.getExistingDirectory(
                self.machine().parent(),
                self.machine().parent().title + ' - Export',
                QtCore.QDir.currentPath())
        else:
            (fileName, _) = QtWidgets.QFileDialog.getSaveFileName(
                self.machine().parent(),
                self.machine().parent().title + ' - Export',
                QtCore.QDir.currentPath() + '/output',
                self.states.edit.infer_dialog_filter())
        if len(fileName) > 0:
            self.data.target = Path(fileName)
            return True
        return False

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
