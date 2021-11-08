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

"""StepAlignOptions, StepAlignSets"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui

from tempfile import TemporaryDirectory
from pathlib import Path

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.threading
import itaxotools.common.resources # noqa

from itaxotools.common.utility import AttrDict

from itaxotools.concatenator import (
    FileType, FileFormat, read_from_path, write_to_path)
from itaxotools.concatenator.library.operators import (
    OpFilterGenes, OpApplyToGene)

from itaxotools.mafftpy import MultipleSequenceAlignment

from .. import model
from .. import widgets
from .. import step_state_machine as ssm

from .wait import StepWaitBar


def work_mafft(input, output, strategy):
    from sys import stdout
    align = MultipleSequenceAlignment(input)
    align.params.general.strategy = strategy
    align.log = stdout
    align.run()
    align.fetch(output)


class RichRadioButton(QtWidgets.QRadioButton):
    def __init__(self, text, desc, parent):
        super().__init__(text, parent)
        self.desc = desc
        self.setStyleSheet("""
            RichRadioButton {
                letter-spacing: 1px;
                font-weight: bold;
            }""")

    def paintEvent(self, event):
        super().paintEvent(event)
        painter = QtGui.QPainter()
        painter.begin(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        font = self.font()
        font.setBold(False)
        font.setLetterSpacing(QtGui.QFont.PercentageSpacing, 0)
        painter.setFont(font)

        width = self.size().width()
        height = self.size().height()
        sofar = super().sizeHint().width()

        rect = QtCore.QRect(sofar, 0, width - sofar, height)
        flags = QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        painter.drawText(rect, flags, self.desc)

        painter.end()

    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        x = event.localPos().x()
        w = self.sizeHint().width()
        if x < w:
            self.setChecked(True)

    def sizeHint(self):
        extra = self.fontMetrics().horizontalAdvance(self.desc)
        size = super().sizeHint()
        size += QtCore.QSize(extra, 0)
        return size


class AlignItem(widgets.ModelItem):
    fields = [
        'display_name',
        'action',
        'samples_len',
        'nucleotides',
        'missing',
        'uniform',
        ]

    def __init__(self, parent, charset: model.Charset):
        super().__init__(parent, charset)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemNeverHasChildren)
        self.samples_len = len(self.model.samples)
        self.tag_new('aligned')
        self.refresh()

    def align(self):
        self.aligned = True
        self.refresh()

    def clear(self):
        self.aligned = False
        self.refresh()

    def toggle(self):
        self.aligned = not self.aligned
        self.refresh()

    def refresh(self):
        self.updateField('action')
        self.tag_set('aligned', self.aligned)
        self.treeWidget().signalTagUpdate.emit()
        self.setBold(self.aligned)

    def setBold(self, value):
        font = self.font(0)
        font.setBold(value)
        self.setFont(0, font)

    @property
    def action(self):
        return 'Align' if self.aligned else '-'


class StepAlignOptions(ssm.StepState):

    title = 'Align Sequences'
    description = 'Select alignment strategy'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = AttrDict()
        self.data.skip = False
        self.data.last = None

    def draw(self):
        widget = QtWidgets.QWidget()

        path = common.resources.get(
            'itaxotools.concatenator_gui', 'docs/mafft.html')
        label = widgets.HtmlLabel(path)

        available = QtWidgets.QLabel('Available strategies:')

        auto = RichRadioButton('Auto', '(depends on data size)', widget)
        desc = '(fast, progressive method, recommended for >2,000 sequences)'
        fftns1 = RichRadioButton('FFT-NS-1', desc, widget)
        desc = '(thorough, slow, recommended for <200 sequences)'
        ginsi = RichRadioButton('G-INS-i', desc, widget)
        skip = QtWidgets.QRadioButton('Skip alignment', widget)
        skip.setStyleSheet("QRadioButton { letter-spacing: 1px; }")
        auto.setChecked(True)

        radios = QtWidgets.QVBoxLayout()
        radios.addWidget(auto)
        radios.addWidget(fftns1)
        radios.addWidget(ginsi)
        radios.addSpacing(16)
        radios.addWidget(skip)
        radios.setSpacing(16)
        radios.setContentsMargins(16, 0, 0, 0)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(label)
        layout.addSpacing(8)
        layout.addWidget(available)
        layout.addLayout(radios)
        layout.addStretch(1)
        layout.setSpacing(24)
        layout.setContentsMargins(0, 0, 0, 32)
        widget.setLayout(layout)

        self.radios = AttrDict()
        self.radios.auto = auto
        self.radios.fftns1 = fftns1
        self.radios.ginsi = ginsi
        self.skip = skip

        return widget

    def get_strategy(self):
        for strategy in self.radios:
            if self.radios[strategy].isChecked():
                return strategy

    def onExit(self, event):
        super().onExit(event)
        self.data.skip = self.skip.isChecked()
        strategy = self.get_strategy()
        if self.data.last != strategy:
            self.data.last = strategy
            self.timestamp_set()


class StepAlignSetsEdit(ssm.StepTriStateEdit):

    description = 'Select which character sets to align'

    def onEntry(self, event):
        super().onEntry(event)
        last_filter_update = self.machine().states.filter.timestamp_get()
        if last_filter_update > self.timestamp_get():
            self.populate_view()
            self.timestamp_set()

    def populate_view(self):
        self.view.clear()
        charsets = [
            cs for cs in self.machine().states.input.data.charsets.values()
            if cs.translation is not None]
        for charset in charsets:
            AlignItem(self.view, charset)
        self.view.resizeColumnsToContents()
        self.sets.setValue(len(charsets))
        self.updateSummary()

    def draw(self):
        widget = QtWidgets.QWidget()

        text = (
            'Mark character sets for alignment by selecting them'
            'and then clicking "Align".')
        label = QtWidgets.QLabel(text)

        frame = self.draw_frame()

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(label)
        layout.addWidget(frame, 1)
        layout.setSpacing(24)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)

        return widget

    def draw_summary(self):
        sets = widgets.InfoLabel('Total Sets')
        marked = widgets.InfoLabel('Marked', 0)

        sets.setToolTip('Total number of character sets.')
        marked.setToolTip('Number of character sets pending alignment.')

        summary = QtWidgets.QHBoxLayout()
        summary.addWidget(sets)
        summary.addWidget(marked)
        summary.addStretch(1)
        summary.setSpacing(24)
        summary.setContentsMargins(4, 0, 4, 0)

        self.sets = sets
        self.marked = marked

        return summary

    def draw_frame(self):
        frame = common.widgets.Frame()

        view = widgets.TreeWidget()
        view.signalTagUpdate.connect(self.updateSummary)
        view.itemActivated.connect(self.handleActivated)
        view.setIndentation(0)
        view.setColumnCount(6, 2)
        view.setHeaderLabels([
            'Name', 'Action', 'Samples',
            'Nucleotides', 'Missing', 'Uniform'])

        headerItem = view.headerItem()
        headerItem.setToolTip(0, 'Character set name')
        headerItem.setToolTip(1, 'Pending action')
        headerItem.setToolTip(2, 'Total number of sequences')
        headerItem.setToolTip(3, 'Total number of nucleotide characters')
        headerItem.setToolTip(4, 'Proportion of missing data')
        headerItem.setToolTip(5, 'Are all sequences of the same length?')

        all = common.widgets.PushButton('Align All', onclick=self.handleAll)
        align = common.widgets.PushButton('Align ', onclick=self.handleAlign)
        clear = common.widgets.PushButton('Clear', onclick=self.handleClear)

        search = widgets.ViewSearchWidget(self, view)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(all)
        controls.addWidget(align)
        controls.addWidget(clear)
        controls.addStretch(1)
        controls.addWidget(search)
        controls.setSpacing(8)
        controls.setContentsMargins(0, 0, 0, 0)

        summary = self.draw_summary()

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(controls)
        layout.addWidget(view, 1)
        layout.addLayout(summary)
        layout.setSpacing(12)
        layout.setContentsMargins(16, 16, 16, 12)
        frame.setLayout(layout)

        self.view = view
        self.frame = frame
        self.search = search

        return frame

    def updateSummary(self):
        self.marked.setValue(self.view.tag_get('aligned'))

    def handleAlign(self, checked=False):
        item = None
        for item in self.view.selectedItems():
            item.align()
        self.view.scrollToItem(item)

    def handleClear(self, checked=False):
        item = None
        for item in self.view.selectedItems():
            item.clear()
        self.view.scrollToItem(item)

    def handleAll(self, checked=False):
        self.view.selectAll()
        self.handleAlign()

    def handleActivated(self, item, column):
        item.toggle()
        self.view.scrollToItem(item)


class StepAlignSetsWait(StepWaitBar):
    pass


class StepAlignSetsDone(ssm.StepTriStateDone):
    description = 'Alignment complete'

    def onEntry(self, event):
        super().onEntry(event)
        if self.result:
            s = 's' if self.result > 1 else ''
            self.parent().update(
                text=f'Successfully aligned {str(self.result)} sequence{s}.')
        else:
            self.parent().update(
                text='Successfully aligned sequences.')


class StepAlignSetsFail(ssm.StepTriStateFail):
    description = 'Alignment failed'

    def onEntry(self, event):
        super().onEntry(event)
        message = f'{type(self.exception).__name__}'
        self.parent().update(
            text=f'Sequence alignment failed: {message}')


class StepAlignSets(ssm.StepTriState):
    title = 'Align Sequences'

    StepEdit = StepAlignSetsEdit
    StepWait = StepAlignSetsWait
    StepDone = StepAlignSetsDone
    StepFail = StepAlignSetsFail

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.temp_prep = None
        self.temp_cache = None
        self.charsets_cached = set()
        self.process = None

    def onEntry(self, event):
        super().onEntry(event)
        input_update = self.machine().states.input.timestamp_get()
        options_update = self.machine().states.align_options.timestamp_get()
        last_update = max(input_update, options_update)
        if last_update > self.timestamp_get():
            self.temp_cache = TemporaryDirectory(prefix='concat_mafft_cache_')
            self.charsets_cached = set()

    def work(self):
        charsets = {
            k for k, v in self.machine().states.input.data.charsets.items()
            if v.aligned and v.translation is not None
            and v.name not in self.charsets_cached}
        files = [
            f for f in self.machine().states.input.data.files.values()
            if any(cs.name in charsets for cs in f.charsets.values())]
        with self.states.wait.redirect():
            self.work_prep(files, charsets)
            self.work_align(charsets)
        self.timestamp_set()
        return len(charsets)

    def work_prep(self, files, charsets):
        def checker_func(series):
            self.worker.check()
            return series

        self.temp_prep = TemporaryDirectory(prefix='concat_mafft_prep_')
        path = Path(self.temp_prep.name)
        print('Preparing files for selected charsets...\n')
        self.update(0, 0, 'Preparing files...')
        for count, file in enumerate(files):
            text = f'Preparing file {count + 1}/{len(files)}: {file.name}'
            print(text)
            self.update(0, 0, text)
            stream = read_from_path(file.path)
            stream = (
                stream
                .pipe(OpApplyToGene(checker_func))
                .pipe(OpFilterGenes(charsets))
                )
            write_to_path(stream, path,FileType.Directory, FileFormat.Fasta)
        self.worker.check()
        print('\nDone preparing files')
        print(f'\n{"-"*20}\n')

    def work_align(self, charsets):
        # We are currently starting double the threads than necessary
        # We could also benefit from the use of a process pool
        # and running the mafft workers in parallel
        strategy = self.machine().states.align_options.get_strategy()
        print(f'Starting alignment for {len(charsets)} sequences '
              f'({len(self.charsets_cached)} already cached)...')
        print(f'Selected strategy: {strategy}\n')
        files = Path(self.temp_prep.name).glob('*')
        outdir = Path(self.temp_cache.name)
        for count, input in enumerate(files):
            charset = input.stem
            text = f'Aligning sequence {count + 1}/{len(charsets)}: {charset}'
            self.update(count, len(charsets), text)
            print(text + '\n')
            output = outdir / input.name

            loop = QtCore.QEventLoop()
            self.process = common.threading.Process(
                work_mafft, input, output, strategy)
            self.process.setStream(self.states.wait.logio)
            self.process.done.connect(loop.quit)
            self.process.fail.connect(self.worker.fail)
            self.process.start()
            loop.exec()

            self.worker.check()
            self.states.wait.reset()
            self.charsets_cached.add(charset)
            print(f'\nAligned {charset}')
            print(f'\n{"-"*20}\n')
        self.update(1, 1, 'Done')
        print('Completed sequence alignment.')

    def onCancel(self, exception):
        self.process.quit()

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

    def skipAll(self):
        skip = self.machine().states['align_options'].data.skip
        return bool(skip)

    def skipWait(self):
        skip = self.states['edit'].marked.value == 0
        return bool(skip)

    def filterCancel(self, event):
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Question)
        msgBox.setText('Cancel alignment?')
        msgBox.setStandardButtons(
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.No)
        res = self.machine().parent().msgShow(msgBox)
        return res == QtWidgets.QMessageBox.Yes
