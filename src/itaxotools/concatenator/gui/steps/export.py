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

from dataclasses import dataclass
from lorem_text import lorem
from random import randint

import enum

from .. import step_state_machine as ssm

from .wait import StepWaitBar


def dummy_work(state, count, max, lines, period):
    print('')
    while True:
        # if count == 101:
        #     raise Exception('ohno')
        now = lorem.words(3)
        print(f'\nFile {count}/{max} {now}')
        for i in range(1, lines):
            print(lorem.words(randint(3, 12)))
        text = f'File {count}/{max}: {now}'
        state.update(count, max, text)
        if count >= max:
            break
        for i in range(0, int(period/10)):
            QtCore.QThread.msleep(10)
            state.worker.check()
        count += 1
    return max


@dataclass
class FileFormat:
    text: str
    many_files: bool
    order_matters: bool


@dataclass
class FileCompression:
    text: str
    extension: str


class CodonOptions(enum.IntEnum):
    Mark = enum.auto()
    Export = enum.auto()


class OrderingOptions(enum.IntEnum):
    SameAsInput = enum.auto()
    Alphabetical = enum.auto()


class DataObject(object):
    pass


class StepExportEdit(ssm.StepTriStateEdit):

    description = 'Configure output'

    def __init__(self, *args, **kwargs):
        self.fileFormats = [
            FileFormat('Interleaved NEXUS', False, True),
            FileFormat('Concatenated Fasta', False, True),
            FileFormat('Concatenated Phylip', True, True),
            FileFormat('Distributed Fasta', False, False),
            FileFormat('Ali Alignment', True, False),
            FileFormat('Partitionfinder', True, True),
            ]
        self.fileCompressions = [
            FileCompression('None', None),
            FileCompression('Zip', 'zip'),
            FileCompression('Tar', 'tar'),
            ]
        self.data = DataObject()
        super().__init__(*args, **kwargs)

    def draw(self):
        widget = QtWidgets.QWidget()

        text = ('Quisque tortor est, porttitor sed viverra ut, '
                'pharetra at nunc. Aenean vel congue dui. '
                'Vivamus auctor, quam se. \n'
                'Lorem ipsum dolor sit amet, consectetur adipiscing elit.'
                )
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

        format = QtWidgets.QComboBox()
        order = QtWidgets.QComboBox()
        compression = QtWidgets.QComboBox()

        format.setToolTip((
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
        label_order = QtWidgets.QLabel('Character set order:')
        label_compression = QtWidgets.QLabel('File compression:')

        label_format.setToolTip(format.toolTip())
        label_order.setToolTip(order.toolTip())
        label_compression.setToolTip(compression.toolTip())


        for fileFormat in self.fileFormats:
            format.addItem(fileFormat.text, fileFormat)
        for fileCompression in self.fileCompressions:
            compression.addItem(fileCompression.text, fileCompression)
        order.addItem('Alphabetical', OrderingOptions.Alphabetical)
        order.addItem('Same as input', OrderingOptions.SameAsInput)

        format.activated.connect(self.formatChanged)

        layout.addWidget(label_format, 0, 0)
        layout.addWidget(label_order, 1, 0)
        layout.addWidget(label_compression, 2, 0)
        layout.addWidget(format, 0, 1)
        layout.addWidget(order, 1, 1)
        layout.addWidget(compression, 2, 1)

        layout.setColumnMinimumWidth(1, 160)
        layout.setSpacing(16)
        layout.setContentsMargins(24, 16, 24, 16)

        self.format = format
        self.order = order
        self.compression = compression

        return layout

    def draw_right(self):
        layout = QtWidgets.QVBoxLayout()

        mark = QtWidgets.QRadioButton(
            'Write codon sets in output information block.')
        mark.setToolTip((
            'Codon sets will be written at:' + '\n'
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

        layout.addWidget(mark)
        layout.addWidget(export)
        layout.addSpacing(16)
        layout.addWidget(exclude)
        layout.setSpacing(8)
        layout.setContentsMargins(24, 16, 24, 16)

        self.codonOptions = options
        self.exclude = exclude

        return layout

    def formatChanged(self, index):
        format = self.format.currentData()
        self.order.setEnabled(format.order_matters)


class StepExportWait(StepWaitBar):
    def onEntry(self, event):
        super().onEntry(event)
        for i in range(1, 10):
            self.logio.writeline(lorem.words(randint(3, 12)))


class StepExportDone(ssm.StepTriStateDone):

    def onEntry(self, event):
        repeat = ssm.NavigateEvent(ssm.NavigateEvent.Event.Next)
        self.machine().postEvent(repeat)


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

    def work(self):
        with self.states['wait'].redirect():
            return dummy_work(self, 2, 10, 2, 80)

    def filterNext(self):
        return True
