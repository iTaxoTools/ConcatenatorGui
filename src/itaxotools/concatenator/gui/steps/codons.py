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

"""StepCodons"""

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


class ItemDelegate(QtWidgets.QStyledItemDelegate):
    def __init__(self, view, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.view = view

    def createEditor(self, parent, option, index):
        if index.column() == 0:
            return super().createEditor(parent, option, index)
        else:
            item = self.view.itemFromIndex(index)
            self.view.itemActivated.emit(item, index.column())
            return None


class CodonItem(widgets.WidgetItem):
    fields = [
        'name',
        'action',
        'frame',
        'code',
        'samples',
        'nucleotides',
        'missing',
        ]
    actions = {
        'Subset': 'subsetted',
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemIsEditable |
                      QtCore.Qt.ItemNeverHasChildren)
        self.file = None
        self.name = lorem.words(randint(2, 6)).replace(' ', '_')
        self.action = '-'
        self.frame = ''
        self.code = ''
        self.samples = randint(6, 300)
        self.nucleotides = randint(200, 3000000)
        self.uniform = ['Yes', 'No'][randint(0, 1)]
        self.missing = randint(0, 9999) / 10000

    def subset(self):
        if self.action == 'Subset':
            return
        self.setAction('Subset')
        self.setBold(True)
        self.frame = 'Auto'
        self.code = 'Auto'

    def clear(self):
        self.setAction('-')
        self.setBold(False)
        self.frame = ''
        self.code = ''

    def toggle(self):
        if self.action == '-':
            self.subset()
        else:
            self.clear()

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


class TreeWidget(widgets.TreeWidget):
    signalSummaryUpdate = QtCore.Signal(str, int)


class DataObject(object):
    pass


class StepCodonsEdit(ssm.StepSubState):

    description = 'Select which character sets to subset'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = DataObject()
        self.add_dummy_contents()

    def add_dummy_contents(self):
        count = randint(500, 2000)
        for i in range(0, count):
            CodonItem(self.view)
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
        subsetted = widgets.InfoLabel('Subsetted', 0)

        for item in [sets, subsetted]:
            item.setToolTip(lorem.words(randint(5, 15)))

        summary = QtWidgets.QHBoxLayout()
        summary.addWidget(sets)
        summary.addWidget(subsetted)
        summary.addStretch(1)
        summary.setSpacing(24)
        summary.setContentsMargins(4, 0, 4, 0)

        self.sets = sets
        self.subsetted = subsetted

        return summary

    def draw_frame(self):
        frame = common.widgets.Frame()

        view = TreeWidget()
        view.setItemDelegate(ItemDelegate(view))
        view.signalSummaryUpdate.connect(self.handleSummaryUpdate)
        view.itemActivated.connect(self.handleActivated)
        view.setIndentation(0)
        view.setColumnCount(7, 4)
        view.setHeaderLabels([
            'Name', 'Action', 'Frame', 'Code', 'Samples',
            'Nucleotides', 'Missing'])

        headerItem = view.headerItem()
        headerItem.setToolTip(0, lorem.words(13))
        for col in range(1, 6):
            headerItem.setToolTip(col, lorem.words(randint(5, 15)))

        subset = common.widgets.PushButton('Subset', onclick=self.handleSubset)
        options = common.widgets.PushButton('Options', onclick=self.handleOptions)
        clear = common.widgets.PushButton('Clear', onclick=self.handleClear)

        search = widgets.ViewSearchWidget(self, view)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(subset)
        controls.addWidget(options)
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

    def handleSubset(self, checked=False):
        for item in self.view.selectedItems():
            item.subset()
        self.view.resizeColumnsToContents()

    def handleClear(self, checked=False):
        for item in self.view.selectedItems():
            item.clear()
        self.view.resizeColumnsToContents()

    def handleOptions(self, checked=False):
        print('handleOptions')

    def handleActivated(self, item, column):
        item.toggle()
        self.view.resizeColumnsToContents()

    def handleSummaryUpdate(self, field, change):
        item = getattr(self, field)
        item.setValue(item.value + change)


class StepCodons(ssm.StepTriState):
    title = 'Subset by Codons'

    StepEdit = StepCodonsEdit

    class StepFail(ssm.StepSubState):
        description = 'Task failed'

    def onFail(self, exception):
        raise exception
