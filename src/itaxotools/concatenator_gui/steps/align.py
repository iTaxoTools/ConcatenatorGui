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

from lorem_text import lorem
from random import randint

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from .. import widgets
from .. import step_state_machine as ssm

from .wait import StepWaitBar


def dummy_work(state, count, max, lines, period):
    print('')
    while True:
        # if count == 101:
        #     raise Exception('ohno')
        now = lorem.words(3)
        print(f'\nStep {count}/{max} {now}')
        for i in range(1, lines):
            print(lorem.words(randint(3, 12)))
        text = f'Sequence {count}/{max}: {now}'
        state.update(count, max, text)
        if count >= max:
            break
        for i in range(0, int(period/10)):
            QtCore.QThread.msleep(10)
            state.worker.check()
        count += 1
    return max


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


class AlignItem(widgets.WidgetItem):
    fields = [
        'name',
        'action',
        'samples',
        'nucleotides',
        'missing',
        'uniform',
        ]
    actions = {
        'Align': 'marked',
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemNeverHasChildren)
        self.file = None
        self.name = lorem.words(randint(2, 6)).replace(' ', '_')
        self.action = '-'
        self.samples = randint(6, 300)
        self.nucleotides = randint(200, 3000000)
        self.uniform = ['Yes', 'No'][randint(0, 1)]
        self.missing = randint(0, 9999) / 10000

    def align(self):
        self.setAction('Align')
        self.setBold(True)

    def clear(self):
        self.setAction('-')
        self.setBold(False)

    def toggle(self):
        if self.action == '-':
            self.align()
        else:
            self.clear()

    def setAction(self, value):
        if self.treeWidget():
            signal = self.treeWidget().signalSummaryUpdate
            if self.action != value:
                if self.action in self.actions:
                    signal.emit(self.actions[self.action], -1)
                if value in self.actions:
                    signal.emit(self.actions[value], 1)
        self.action = value

    def setBold(self, value):
        font = self.font(0)
        font.setBold(value)
        self.setFont(0, font)


class TreeWidget(widgets.TreeWidget):
    signalSummaryUpdate = QtCore.Signal(str, int)


class DataObject(object):
    pass


class StepAlignOptions(ssm.StepState):

    title = 'Align Sequences'
    description = 'Select alignment strategy'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = DataObject()
        self.data.skip = False

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

        self.auto = auto
        self.fftns1 = fftns1
        self.ginsi = ginsi
        self.skip = skip

        return widget

    def onExit(self, event):
        super().onExit(event)
        self.data.skip = self.skip.isChecked()


class StepAlignSetsEdit(ssm.StepTriStateEdit):

    description = 'Select which character sets to align'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = DataObject()
        self.add_dummy_contents()

    def onEntry(self, event):
        super().onEntry(event)
        self.updateFooter()

    def add_dummy_contents(self):
        count = randint(500, 2000)
        for i in range(0, count):
            AlignItem(self.view)
        self.view.resizeColumnsToContents()
        self.sets.setValue(count)

    def draw(self):
        widget = QtWidgets.QWidget()

        text = (
            'Select character sets and click "Align" to mark them '
            'for alignment. When ready, click "Start".')
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

        view = TreeWidget()
        view.signalSummaryUpdate.connect(self.handleSummaryUpdate)
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

    def handleAlign(self, checked=False):
        for item in self.view.selectedItems():
            item.align()
        self.view.scrollToItem(item)

    def handleClear(self, checked=False):
        for item in self.view.selectedItems():
            item.clear()
        self.view.scrollToItem(item)

    def handleAll(self, checked=False):
        self.view.selectAll()
        self.handleAlign()

    def handleActivated(self, item, column):
        item.toggle()
        self.view.scrollToItem(item)

    def handleSummaryUpdate(self, field, change):
        item = getattr(self, field)
        item.setValue(item.value + change)
        if field == 'marked':
            self.updateFooter()

    def updateFooter(self):
        if self.marked.value == 0:
            self.footer.next.setText('&Skip >')
        else:
            self.footer.next.setText('&Start')


class StepAlignSetsWait(StepWaitBar):
    def onEntry(self, event):
        super().onEntry(event)
        for i in range(1, 50):
            self.logio.writeline(lorem.words(randint(3, 12)))


class StepAlignSetsDone(ssm.StepTriStateDone):
    description = 'Alignment complete'

    def onEntry(self, event):
        super().onEntry(event)
        self.parent().update(
            text=f'Successfully aligned {str(self.result)} sequences.')


class StepAlignSetsFail(ssm.StepTriStateFail):
    description = 'Alignment failed'

    def onEntry(self, event):
        super().onEntry(event)
        self.parent().update(
            text=f'Sequence alignment failed: {str(self.exception)}')


class StepAlignSets(ssm.StepTriState):
    title = 'Align Sequences'

    StepEdit = StepAlignSetsEdit
    StepWait = StepAlignSetsWait
    StepDone = StepAlignSetsDone
    StepFail = StepAlignSetsFail

    def work(self):
        with self.states['wait'].redirect():
            return dummy_work(self, 42, 200, 10, 20)

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
        res = msgBox.exec()
        return res == QtWidgets.QMessageBox.Yes
