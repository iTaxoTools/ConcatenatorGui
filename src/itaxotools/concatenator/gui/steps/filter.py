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

"""StepFilter"""

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


class ItemDelegate(QtWidgets.QStyledItemDelegate):
    def createEditor(self, parent, option, index):
        if index.column() == 0:
            return super().createEditor(parent, option, index)
        else:
            return None


class FilterItem(widgets.WidgetItem):
    fields = [
        'name',
        'action',
        'samples',
        'nucleotides',
        'uniform',
        'missing',
        ]
    actions = {
        'Delete': 'deleted',
        'Rename': 'renamed',
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemIsEditable |
                      QtCore.Qt.ItemNeverHasChildren)
        self.file = None
        alignment = QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter
        for col in range(2, 6):
            self.setTextAlignment(col, alignment)
        self.name = lorem.words(randint(2, 6)).replace(' ', '_')
        self.name_original = self.name
        self.action = '-'
        self.samples = randint(6, 300)
        self.nucleotides = randint(200, 3000000)
        self.uniform = ['Yes', 'No'][randint(0, 1)]
        self.missing = randint(0, 9999) / 10000

    def setData(self, column, role, value):
        if not value:
            return
        super().setData(column, role, value)
        if role != QtCore.Qt.ItemDataRole.EditRole:
            return
        if column > 0:
            raise RuntimeError(f'Cannot edit FilterItem column: {column}')
        if value != self.name:
            self.rename(value)

    def rename(self, value):
        self.name = value
        if self.name != self.name_original:
            self.setAction('Rename')
            self.setBold(True)
        else:
            self.setAction('-')
            self.setBold(False)

    def delete(self):
        self.setAction('Delete')
        self.setBold(True)

    def clear(self):
        self.name = self.name_original
        self.setAction('-')
        self.setBold(False)

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


class StepFilter(ssm.StepState):

    title = 'Filter Character Sets'
    description = 'Rename or delete character sets'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = DataObject()
        self.add_dummy_contents()

    def add_dummy_contents(self):
        count = randint(500, 2000)
        for i in range(0, count):
            FilterItem(self.view)
        self.view.resizeColumnsToContents()
        self.sets.setValue(count)

    def work(self):
        time = randint(500, 2000)
        for i in range(0, int(time/10)):
            QtCore.QThread.msleep(10)
            self.worker.check()

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
        renamed = widgets.InfoLabel('Renamed', 0)
        deleted = widgets.InfoLabel('Deleted', 0)

        for item in [sets, renamed, deleted]:
            item.setToolTip(lorem.words(randint(5, 15)))

        summary = QtWidgets.QHBoxLayout()
        summary.addWidget(sets)
        summary.addWidget(renamed)
        summary.addWidget(deleted)
        summary.addStretch(1)
        summary.setSpacing(24)
        summary.setContentsMargins(4, 0, 4, 0)

        self.sets = sets
        self.renamed = renamed
        self.deleted = deleted

        return summary

    def draw_frame(self):
        frame = common.widgets.Frame()

        view = TreeWidget()
        view.setItemDelegate(ItemDelegate())
        view.signalSummaryUpdate.connect(self.handleSummaryUpdate)
        view.itemActivated.connect(self.handleActivated)
        view.itemDoubleClicked.connect(self.handleActivated)
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

        rename = common.widgets.PushButton('Rename', onclick=self.handleRename)
        delete = common.widgets.PushButton('Delete', onclick=self.handleDelete)
        clear = common.widgets.PushButton('Clear', onclick=self.handleClear)

        QtGui.QShortcut(QtGui.QKeySequence('F2'), view, self.handleRename)
        QtGui.QShortcut(QtGui.QKeySequence.Delete, view, self.handleDelete)
        QtGui.QShortcut(QtGui.QKeySequence.Undo, view, self.handleClear)

        search = widgets.ViewSearchWidget(self, view)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(rename)
        controls.addWidget(delete)
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
        self.delete = delete
        self.search = search

        return frame

    def handleRename(self, checked=False):
        index = self.view.currentIndex().siblingAtColumn(0)
        self.view.edit(index)
        self.view.resizeColumnToContents(1)

    def handleDelete(self, checked=False):
        for item in self.view.selectedItems():
            item.delete()
        self.view.resizeColumnToContents(1)

    def handleClear(self, checked=False):
        for item in self.view.selectedItems():
            item.clear()
        self.view.resizeColumnToContents(1)

    def handleActivated(self, item, column):
        index = self.view.indexFromItem(item).siblingAtColumn(0)
        self.view.edit(index)

    def handleSummaryUpdate(self, field, change):
        item = getattr(self, field)
        item.setValue(item.value + change)
