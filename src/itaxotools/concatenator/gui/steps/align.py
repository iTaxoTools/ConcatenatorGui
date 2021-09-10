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

import sys

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from .. import widgets
from .. import step_state_machine as ssm


def dummy_work(self, count, max, lines, period):
    with common.io.redirect(sys, 'stdout', self.logio):
        print('')
        while True:
            text = lorem.words(3)
            print(f'\nStep {count}/{max} {text}')
            for i in range(1, lines):
                print(lorem.words(randint(3, 12)))
            self.updateSignal.emit((count, max, text))
            if count >= max:
                break
            for i in range(0, int(period/10)):
                QtCore.QThread.msleep(10)
                self.worker.check()
            count += 1


class AlignItem(widgets.WidgetItem):
    fields = [
        'name',
        'action',
        'samples',
        'nucleotides',
        'uniform',
        'missing',
        ]
    actions = {
        'Align': 'aligned',
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemNeverHasChildren)
        self.file = None
        alignment = QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter
        for col in range(2, 6):
            self.setTextAlignment(col, alignment)
        self.name = lorem.words(randint(2, 6)).replace(' ', '_')
        self.action = '-'
        self.samples = randint(6, 300)
        self.nucleotides = randint(200, 3000000)
        self.uniform = ['Yes', 'No'][randint(0, 1)]
        self.missing = randint(0, 9999) / 10000

    def align(self):
        self.setAction('Align')
        self.setBold(True)
        self.treeWidget().resizeColumnToContents(1)

    def ignore(self):
        self.setAction('-')
        self.setBold(False)
        self.treeWidget().resizeColumnToContents(1)

    def toggle(self):
        if self.action == '-':
            self.align()
        else:
            self.ignore()

    def setAction(self, value):
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
        self.setFont(1, font)


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

        text = ('Quisque tortor est, porttitor sed viverra ut, '
                'pharetra at nunc. Aenean vel congue dui. '
                'Vivamus auctor, quam se. \n'
                'Lorem ipsum dolor sit amet, consectetur adipiscing elit.'
                )
        label = QtWidgets.QLabel(text)

        placeholder = QtWidgets.QLabel('Strategy options')

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(label)
        layout.addWidget(placeholder)
        layout.addStretch(1)
        layout.setSpacing(24)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)

        return widget


class StepAlignSetsEdit(ssm.StepSubState):

    description = 'Select which character sets to align'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = DataObject()
        self.add_dummy_contents()

    def add_dummy_contents(self):
        count = randint(5, 50)
        for i in range(0, count):
            AlignItem(self.view)
        self.view.resizeColumnsToContents()
        self.sets.setValue(count)

    def draw(self):
        widget = QtWidgets.QWidget()

        text = ('Quisque tortor est, porttitor sed viverra ut, '
                'pharetra at nunc. Aenean vel congue dui. '
                'Vivamus auctor, quam se. \n'
                'Lorem ipsum dolor sit amet, consectetur adipiscing elit.'
                )
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
        aligned = widgets.InfoLabel('Aligned', 0)

        for item in [sets, aligned]:
            item.setToolTip(lorem.words(randint(5, 15)))

        summary = QtWidgets.QHBoxLayout()
        summary.addWidget(sets)
        summary.addWidget(aligned)
        summary.addStretch(1)
        summary.setSpacing(24)
        summary.setContentsMargins(4, 0, 4, 0)

        self.sets = sets
        self.aligned = aligned

        return summary

    def draw_frame(self):
        frame = common.widgets.Frame()

        view = TreeWidget()
        view.signalSummaryUpdate.connect(self.handleSummaryUpdate)
        view.itemActivated.connect(self.handleActivated)
        view.setIndentation(0)
        view.setColumnCount(6)
        view.setHeaderLabels([
            ' Name', ' Action', ' Samples',
            ' Nucleotides', ' Uniform', ' Missing'])

        headerItem = view.headerItem()
        headerItem.setToolTip(0, lorem.words(13))
        for col in range(1, 6):
            headerItem.setToolTip(col, lorem.words(randint(5, 15)))
        alignment = QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        headerItem.setTextAlignment(1, alignment)

        align = common.widgets.PushButton('Align', onclick=self.handleAlign)
        all = common.widgets.PushButton('Align All', onclick=self.handleAll)
        ignore = common.widgets.PushButton('Ignore', onclick=self.handleIgnore)

        search = widgets.ViewSearchWidget(self, view)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(align)
        controls.addWidget(all)
        controls.addWidget(ignore)
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
        self.frame = align
        self.delete = ignore
        self.search = all

        return frame

    def handleAlign(self, checked=False):
        for item in self.view.selectedItems():
            item.align()

    def handleIgnore(self, checked=False):
        for item in self.view.selectedItems():
            item.ignore()

    def handleAll(self, checked=False):
        self.view.selectAll()
        self.handleAlign()

    def handleActivated(self, item, column):
        item.toggle()

    def handleSummaryUpdate(self, field, change):
        item = getattr(self, field)
        item.setValue(item.value + change)


class StepAlignSetsWait(ssm.StepSubState):
    description = 'Please wait...'

    def draw(self):
        widget = QtWidgets.QWidget()

        progtext = QtWidgets.QLabel('No task')
        progtext.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        progbar = QtWidgets.QProgressBar()
        progbar.setTextVisible(False)
        progbar.setStyleSheet("""
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
        logger = common.widgets.TextEditLogger()
        logger.setFont(QtGui.QFontDatabase.systemFont(
            QtGui.QFontDatabase.FixedFont))
        logger.document().setDocumentMargin(10)
        logio = common.io.TextEditLoggerIO(logger)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(progtext)
        layout.addSpacing(6)
        layout.addWidget(progbar)
        layout.addSpacing(16)
        layout.addWidget(logger, 1)
        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)

        self.parent().logio = logio
        self.parent().logger = logger
        self.parent().progtext = progtext
        self.parent().progbar = progbar

        return widget

    def onEntry(self, event):
        super().onEntry(event)
        self.parent().logger.clear()
        for i in range(1, 50):
            self.parent().logio.writeline(lorem.words(randint(3, 12)))


class StepAlignSets(ssm.StepTriState):
    title = 'Align Sequences'

    StepEdit = StepAlignSetsEdit
    StepWait = StepAlignSetsWait

    class StepFail(ssm.StepSubState):
        description = 'Task failed'

    def onFail(self, exception):
        raise exception

    def onEntry(self, event):
        skip = self.machine().states['align_options'].data.skip
        if skip:
            repeat = ssm.NavigateEvent(event.event)
            self.machine().postEvent(repeat)
        else:
            super().onEntry(event)

    def updateProgress(self, tuple):
        x, m, w = tuple
        self.progtext.setText(f'Sequence {x}/{m}: {w}')
        self.progbar.setMaximum(m)
        self.progbar.setValue(x)

    def work(self):
        skip = self.machine().states['align_options'].data.skip
        skip = skip or self.states['edit'].aligned.value == 0
        if not skip:
            dummy_work(self, 42, 200, 10, 20)

    def filterCancel(self):
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle('Cancel')
        msgBox.setIcon(QtWidgets.QMessageBox.Question)
        msgBox.setText('Cancel current task?')
        msgBox.setStandardButtons(
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.No)
        res = msgBox.exec()
        return res == QtWidgets.QMessageBox.Yes
