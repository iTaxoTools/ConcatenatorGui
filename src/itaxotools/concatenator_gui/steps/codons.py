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

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from .. import widgets
from .. import step_state_machine as ssm

from .wait import StepWaitBar


def dummy_work(state, count, max, lines, period):
    print('')
    while True:
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


class CodonItem(widgets.WidgetItem):
    fields = [
        'name',
        'action',
        'code',
        'frame',
        'naming',
        ]
    actions = {
        'Split': 'marked',
        }
    values = {
        'naming': ['Default', 'Custom'],
        'frame': {
            'Auto-detect': 'Auto-detect',
            '+1': '+1 (starts with 1st codon position)',
            '+2': '+2 (starts with 2nd codon position)',
            '+3': '+3 (starts with 3rd codon position)',
            '-1': '-1 (reverse complement of +1)',
            '-2': '-2 (reverse complement of +2)',
            '-3': '-3 (reverse complement of +3)',
            },
        'code': [
            'Auto-detect',
            'Standard',
            'Vertebrate Mitochondrial',
            'Yeast Mitochondrial',
            'Invertebrate Mitochondrial',
            ],
        }
    defaults = {
        1: '**_1st',
        2: '**_2nd',
        3: '**_3rd',
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemIsEditable |
                      QtCore.Qt.ItemNeverHasChildren)
        self.file = None
        self.name = lorem.words(randint(2, 6)).replace(' ', '_')
        self.action = None
        self.clear()

    def subset(self):
        if self.action == 'Split':
            return
        self.setAction('Split')
        self.setBold(True)
        self.naming = 'Default'
        self.frame = 'Auto-detect'
        self.code = 'Auto-detect'

    def clear(self):
        self.setAction('-')
        self.setBold(False)
        self.naming = ''
        self.frame = ''
        self.code = ''
        self.names = dict(self.defaults)

    def toggle(self):
        if self.action == '-':
            self.subset()
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

    def setData(self, column, role, value):
        super().setData(column, role, value)
        if role == QtCore.Qt.EditRole and column == self.map['naming']:
            if value == 'Default':
                self.setNames()

    def setNames(self, names={1: None, 2: None, 3: None}):
        for key, name in names.items():
            if name:
                self.names[key] = name
            else:
                self.names[key] = self.defaults[key]
        matches = [name for key, name in self.names.items()
                   if name == self.defaults[key]]
        if len(matches) == len(self.defaults):
            self.naming = 'Default'
        else:
            self.naming = 'Custom'

    @classmethod
    def populateComboBox(cls, combo, field):
        values = cls.values[field]
        if isinstance(values, list):
            for value in values:
                combo.addItem(value, value)
        elif isinstance(values, dict):
            for key, value in values.items():
                combo.addItem(value, key)


class PopupItemDelegate(QtWidgets.QStyledItemDelegate):
    def sizeHint(self, option, index):
        return self.parent().size()


class EditComboBox(QtWidgets.QComboBox):
    def __init__(self, delegate, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.delegate = delegate
        self.activated.connect(self.commit)
        self.setView(QtWidgets.QListView(self))
        self.view().setItemDelegate(PopupItemDelegate(self))
        self.view().setStyleSheet("""
        QListView::item {
            background-color: Palette(Light);
            color: Palette(Text);
            padding: 6px;
        }
        QListView::item:selected {
            background-color: Palette(Highlight);
            color: Palette(Light);
        }
        """)

    def event(self, event):
        if isinstance(event, QtGui.QKeyEvent):
            if event.matches(QtGui.QKeySequence.Cancel):
                self.delegate.closeEditor.emit(self)
        return super().event(event)

    def commit(self):
        self.delegate.commitData.emit(self)
        self.delegate.closeEditor.emit(self)

    def paintEvent(self, event):
        # Never draw the combobox itself
        pass


class CodonItemDelegate(QtWidgets.QStyledItemDelegate):
    def __init__(self, view, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.view = view

    def updateEditorGeometry(self, editor, option, index):
        # Center the combobox popup. Ugly but it works.
        # For QTreeView::item { padding: 2px 6px; }
        option.rect += QtCore.QMargins(7, 3, 5, 1)
        super().updateEditorGeometry(editor, option, index)

    def createEditor(self, parent, option, index):
        item = self.view.itemFromIndex(index)
        column = index.column()
        field = item.unmap[column]
        if item.action == '-':
            self.view.itemActivated.emit(item, column)
            return None
        if field in item.values.keys():
            combo = EditComboBox(self, parent)
            item.populateComboBox(combo, field)
            return combo
        else:
            self.view.itemActivated.emit(item, column)
            return None

    def setEditorData(self, editor, index):
        value = index.data(QtCore.Qt.EditRole)
        id = editor.findData(value)
        if (id >= 0):
            editor.setCurrentIndex(id)
        editor.showPopup()

    def setModelData(self, editor, model, index):
        if isinstance(editor, EditComboBox):
            model.setData(index, editor.currentData(), QtCore.Qt.EditRole)
        else:
            super().setModelData(editor, model, index)


class TreeWidget(widgets.TreeWidget):
    signalSummaryUpdate = QtCore.Signal(str, int)


class DataObject(object):
    pass


class OptionsDialog(QtWidgets.QDialog):
    def __init__(self, view, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.view = view

        self.setWindowTitle(self.parent().title)
        self.draw()
        self.setFixedSize(self.sizeHint())

        items = self.view.selectedItems()
        if len(items) == 1:
            item = items[0]
            self.label_char.setText('Character Set:')
            self.char.setText(item.name)
        elif len(items) > 1:
            self.label_char.setText('Character Sets:')
            self.char.setText(f'{len(items)} items selected')

            item = CodonItem(None)
            item.code = self.get_common_attr(items, 'code')
            item.frame = self.get_common_attr(items, 'frame')
            names = {}
            names[1] = self.get_common_attr(
                items, 'names', func=lambda x, a: getattr(x, a)[1])
            names[2] = self.get_common_attr(
                items, 'names', func=lambda x, a: getattr(x, a)[2])
            names[3] = self.get_common_attr(
                items, 'names', func=lambda x, a: getattr(x, a)[3])
            item.setNames(names)

        id = self.code.findData(item.code)
        if (id >= 0):
            self.code.setCurrentIndex(id)

        id = self.frame.findData(item.frame)
        if (id >= 0):
            self.frame.setCurrentIndex(id)

        self.first.setText(item.names[1])
        self.second.setText(item.names[2])
        self.third.setText(item.names[3])

        for widget in [self.char, self.first, self.second, self.third]:
            widget.setCursorPosition(0)

    def draw(self):
        header = QtWidgets.QLabel(
            'Set the codon subset options for all selected character sets.')
        name_desc = QtWidgets.QLabel(
            'Double asterisks (**) are replaced by the character set name.')
        self.label_char = QtWidgets.QLabel('Character Set:')
        label_code = QtWidgets.QLabel('Genetic Code:')
        label_frame = QtWidgets.QLabel('Reading Frame:')
        label_1 = QtWidgets.QLabel('1st Codon:')
        label_2 = QtWidgets.QLabel('2nd Codon:')
        label_3 = QtWidgets.QLabel('3rd Codon:')

        self.char = QtWidgets.QLineEdit('Lorep_impsum_nocte')
        self.char.setReadOnly(True)
        self.char.setStyleSheet("""
            QLineEdit {
                border: 1px solid Palette(Mid);
                border-radius: 2px;
                padding: 2px
            }""")
        self.first = QtWidgets.QLineEdit('**_1st')
        self.second = QtWidgets.QLineEdit('**_2nd')
        self.third = QtWidgets.QLineEdit('**_3rd')

        self.code = QtWidgets.QComboBox()
        CodonItem.populateComboBox(self.code, 'code')
        self.frame = QtWidgets.QComboBox()
        CodonItem.populateComboBox(self.frame, 'frame')

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

        layout = QtWidgets.QGridLayout()
        layout.addWidget(header, 0, 0, 1, 2)

        layout.setRowMinimumHeight(10, 8)
        layout.addWidget(self.label_char, 11, 0)
        layout.addWidget(self.char, 11, 1)
        layout.addWidget(label_code, 12, 0)
        layout.addWidget(self.code, 12, 1)
        layout.addWidget(label_frame, 13, 0)
        layout.addWidget(self.frame, 13, 1)

        layout.setRowMinimumHeight(20, 16)
        layout.addWidget(label_1, 22, 0)
        layout.addWidget(self.first, 22, 1)
        layout.addWidget(label_2, 23, 0)
        layout.addWidget(self.second, 23, 1)
        layout.addWidget(label_3, 24, 0)
        layout.addWidget(self.third, 24, 1)
        layout.addWidget(name_desc, 29, 0, 1, 2)

        layout.setRowMinimumHeight(30, 24)
        layout.addLayout(buttons, 31, 0, 1, 2)

        layout.setContentsMargins(24, 16, 24, 16)
        self.setLayout(layout)

    def get_common_attr(self, items, attr, func=lambda x, a: getattr(x, a)):
        if not items:
            return None
        res = func(items[0], attr)
        for item in items:
            if func(item, attr) != res:
                return None
        return res

    def accept(self):
        super().accept()
        items = self.view.selectedItems()
        for item in items:
            item.subset()
            item.code = self.code.currentData()
            item.frame = self.frame.currentData()
            names = {
                1: self.first.text(),
                2: self.second.text(),
                3: self.third.text(),
                }
            item.setNames(names)
        self.view.scrollToItem(item)


class StepCodonsEdit(ssm.StepTriStateEdit):

    description = 'Select which character sets to split'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = DataObject()
        self.add_dummy_contents()

    def onEntry(self, event):
        super().onEntry(event)
        self.updateFooter()

    def add_dummy_contents(self):
        count = randint(5, 20)
        for i in range(0, count):
            CodonItem(self.view)
        self.view.resizeColumnsToContents()
        self.sets.setValue(count)

    def draw(self):
        widget = QtWidgets.QWidget()

        text = (
            'Double-click an option field to change it. '
            'Click "Edit" to set codon names and bulk editing.'
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
        marked = widgets.InfoLabel('Marked', 0)

        sets.setToolTip('Total number of character sets.')
        marked.setToolTip('Character sets pending a split.')

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
        view.setItemDelegate(CodonItemDelegate(view))
        view.signalSummaryUpdate.connect(self.handleSummaryUpdate)
        view.itemActivated.connect(self.handleActivated)
        view.setIndentation(0)
        view.setColumnCount(5, 5)
        view.setHeaderLabels([
            'Name', 'Action', 'Genetic Code', 'Reading Frame', 'Naming'])

        headerItem = view.headerItem()
        headerItem.setToolTip(0, 'Character set name')
        headerItem.setToolTip(1, 'Pending action')
        headerItem.setToolTip(2, 'Genetic code variant (translation table)')
        headerItem.setToolTip(3, 'Reading frame')
        headerItem.setToolTip(4, 'Codon naming method')

        subset = common.widgets.PushButton('Subset', onclick=self.handleSubset)
        edit = common.widgets.PushButton('Edit', onclick=self.handleOptions)
        clear = common.widgets.PushButton('Clear', onclick=self.handleClear)

        search = widgets.ViewSearchWidget(self, view)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(subset)
        controls.addWidget(edit)
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
        item = None
        for item in self.view.selectedItems():
            item.subset()
        self.view.scrollToItem(item)

    def handleClear(self, checked=False):
        item = None
        for item in self.view.selectedItems():
            item.clear()
        self.view.scrollToItem(item)

    def handleOptions(self, checked=False):
        items = self.view.selectedItems()
        if not items:
            return
        self.dialog = OptionsDialog(self.view, self.machine().parent())
        self.dialog.setModal(True)
        self.dialog.show()

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


class StepCodonsWait(StepWaitBar):
    def onEntry(self, event):
        super().onEntry(event)
        for i in range(1, 50):
            self.logio.writeline(lorem.words(randint(3, 12)))


class StepCodonsDone(ssm.StepTriStateDone):
    description = 'Codon subsetting complete'

    def onEntry(self, event):
        super().onEntry(event)
        marked = self.parent().states['edit'].marked.value
        s = 's' if marked > 1 else ''
        self.parent().update(
            text=f'Successfully split {marked} character set{s}.')


class StepCodonsFail(ssm.StepTriStateFail):
    description = 'Codon subsetting failed'

    def onEntry(self, event):
        super().onEntry(event)
        self.parent().update(
            text=f'Codon subsetting failed: {str(self.exception)}')


class StepCodons(ssm.StepTriState):
    title = 'Subset by Codons'

    StepEdit = StepCodonsEdit
    StepWait = StepCodonsWait
    StepDone = StepCodonsDone
    StepFail = StepCodonsFail

    class StepFail(ssm.StepSubState):
        description = 'Task failed'

    def onFail(self, exception):
        raise exception

    def work(self):
        with self.states['wait'].redirect():
            return dummy_work(self, 42, 100, 10, 20)

    def skipWait(self):
        skip = self.states['edit'].marked.value == 0
        return bool(skip)

    def filterCancel(self, event):
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Question)
        msgBox.setText('Cancel subsetting?')
        msgBox.setStandardButtons(
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.No)
        res = msgBox.exec()
        return res == QtWidgets.QMessageBox.Yes
