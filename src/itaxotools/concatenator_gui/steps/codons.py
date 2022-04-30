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

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from itaxotools.concatenator.library.model import GeneSeries
from itaxotools.concatenator.library.operators import OpUpdateMetadata
from itaxotools.concatenator.library.codons import (
    GeneticCode, ReadingFrame, _GC_DESCRIPTIONS)

from .. import model
from .. import widgets
from .. import step_state_machine as ssm

from .wait import StepWaitBar


GENETIC_CODES = {0: 'Unknown'}
GENETIC_CODES.update({
    k: f'{GeneticCode(k).name}: {_GC_DESCRIPTIONS[k].name}'
    for k in _GC_DESCRIPTIONS.keys()
})

# READING_FRAMES = {0: 'Auto-detect'}
READING_FRAMES = {}
READING_FRAMES.update({
    frame: str(frame) for frame in ReadingFrame
    if frame != 0
})

DEFAULT_NAMES = {
    i+1: name for i, name in enumerate(GeneSeries.defaults['codon_names'])}


class CodonItem(widgets.ModelItem):
    fields = [
        'display_name',
        'action',
        'code_display',
        'frame_display',
        'naming',
        ]
    values = {
        'naming': ['Default', 'Custom'],
        'frame_display': READING_FRAMES,
        'code_display': GENETIC_CODES
        }
    defaults = DEFAULT_NAMES

    def __init__(self, parent, charset: model.Charset):
        super().__init__(parent, charset)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemIsEditable |
                      QtCore.Qt.ItemNeverHasChildren)
        self.tag_new('split')
        self.clear()
        self.refresh()

    def subset(self):
        if self.split:
            return
        self.split = True
        self.naming = 'Default'
        self.frame = ReadingFrame(1)
        self.code = 0
        self.refresh()

    def clear(self):
        self.split = False
        self.naming = ''
        self.frame = ReadingFrame(1)
        self.code = 0
        self.names = dict(self.defaults)
        self.refresh()

    def toggle(self):
        if not self.split:
            self.subset()
        else:
            self.clear()

    def refresh(self):
        self.updateField('action')
        self.updateField('code_display')
        self.updateField('frame_display')
        self.tag_set('split', self.split)
        if self.treeWidget():
            self.treeWidget().signalTagUpdate.emit()
        self.setBold(self.split)

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

    @property
    def action(self):
        return 'Subset' if self.split else '-'

    @property
    def code_display(self):
        if not self.split:
            return ''
        return GeneticCode(self.code).name

    @code_display.setter
    def code_display(self, value):
        self.code = value

    @property
    def frame_display(self):
        if not self.split:
            return ''
        if not self.frame:
            return 'Auto-detect'
        return ReadingFrame(self.frame).label

    @frame_display.setter
    def frame_display(self, value):
        self.frame = value

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
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.view = view

        self.setWindowTitle(self.parent().title)
        self.draw()
        self.setFixedSize(self.sizeHint())

        items = self.view.selectedItems()
        if len(items) == 1:
            item = items[0]
            self.label_char.setText('Marker Name:')
            self.char.setText(item.name)
        elif len(items) > 1:
            self.label_char.setText('Marker Name:')
            self.char.setText(f'{len(items)} items selected')

            item = CodonItem(None, model.Charset(''))
            item.split = self.get_common_attr(items, 'split')
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

        self.split.setChecked(bool(item.split))

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
        self.label_char = QtWidgets.QLabel('Gene Name:')
        label_code = QtWidgets.QLabel('Genetic Code:')
        label_frame = QtWidgets.QLabel('Reading Frame:')
        label_1 = QtWidgets.QLabel('1st Codon:')
        label_2 = QtWidgets.QLabel('2nd Codon:')
        label_3 = QtWidgets.QLabel('3rd Codon:')
        name_desc = QtWidgets.QLabel(
            'Double asterisks (**) are replaced by the marker name.')

        self.split = QtWidgets.QCheckBox(
            ' Subset codon positions for all selected markers.')
        self.char = QtWidgets.QLineEdit('')
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
        self.code.setMaximumWidth(240)
        CodonItem.populateComboBox(self.code, 'code_display')
        self.frame = QtWidgets.QComboBox()
        CodonItem.populateComboBox(self.frame, 'frame_display')

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
        options.addWidget(self.label_char, 11, 0)
        options.addWidget(self.char, 11, 1)
        options.addWidget(label_code, 12, 0)
        options.addWidget(self.code, 12, 1)
        options.addWidget(label_frame, 13, 0)
        options.addWidget(self.frame, 13, 1)

        options.setRowMinimumHeight(20, 16)
        options.addWidget(label_1, 22, 0)
        options.addWidget(self.first, 22, 1)
        options.addWidget(label_2, 23, 0)
        options.addWidget(self.second, 23, 1)
        options.addWidget(label_3, 24, 0)
        options.addWidget(self.third, 24, 1)
        options.addWidget(name_desc, 29, 0, 1, 2)

        options.setContentsMargins(0, 0, 0, 0)
        options.setRowMinimumHeight(30, 24)

        self.options = QtWidgets.QWidget()
        self.options.setLayout(options)

        self.split.stateChanged.connect(self.options.setEnabled)
        self.split.setChecked(True)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.split)
        layout.addWidget(self.options)
        layout.addLayout(buttons)
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
        if self.split.isChecked():
            for item in items:
                item.subset()
                item.code_display = self.code.currentData()
                item.frame_display = self.frame.currentData()
                names = {
                    1: self.first.text(),
                    2: self.second.text(),
                    3: self.third.text(),
                    }
                item.setNames(names)
        else:
            for item in items:
                item.clear()
        self.view.scrollToItem(item)


class StepCodonsEdit(ssm.StepTriStateEdit):

    description = 'Select which genes to subset'

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
            CodonItem(self.view, charset)
        self.view.resizeColumnsToContents()
        self.sets.setValue(len(charsets))
        self.updateSummary()

    def draw(self):
        widget = QtWidgets.QWidget()

        text = (
            'Select markers for codon subsetting in the character set '
            'specifications of the output file. '
            'Double-click an option field to edit it.')
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
        sets = widgets.InfoLabel('Markers')
        marked = widgets.InfoLabel('Selected', 0)

        sets.setToolTip('Total number of markers.')
        marked.setToolTip('Number of markers selected for subsetting.')

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
        view.signalTagUpdate.connect(self.updateSummary)
        view.itemActivated.connect(self.handleActivated)
        view.setIndentation(0)
        view.setColumnCount(5, 5)
        view.setHeaderLabels([
            'Name', 'Action', 'Genetic Code', 'Reading Frame', 'Naming'])

        headerItem = view.headerItem()
        headerItem.setToolTip(0, 'Marker name')
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

    def updateSummary(self):
        self.marked.setValue(self.view.tag_get('split'))

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


class StepCodonsWait(StepWaitBar):
    pass


class StepCodonsDone(ssm.StepTriStateDone):
    description = 'Codon subsetting complete'

    def onEntry(self, event):
        super().onEntry(event)
        marked = self.parent().states.edit.marked.value
        s = 's' if marked > 1 else ''
        self.parent().update(
            text=f'Successfully subsetted {marked} gene{s}.')


class StepCodonsFail(ssm.StepTriStateFail):
    description = 'Codon subsetting failed'

    def onEntry(self, event):
        super().onEntry(event)
        message = f'{type(self.exception).__name__}'
        self.parent().update(
            text=f'Codon subsetting failed: {message}')


class StepCodons(ssm.StepTriState):
    title = 'Character Sets by Codon'

    StepEdit = StepCodonsEdit
    StepWait = StepCodonsWait
    StepDone = StepCodonsDone
    StepFail = StepCodonsFail

    class StepFail(ssm.StepSubState):
        description = 'Task failed'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.metas = dict()

    def work(self):
        with self.states['wait'].redirect():
            print('Codon subset automation is not supported at this time.')
            print('All reading frames will be exported as is.')
            self.metas = {
                item.name: {
                    'genetic_code': GeneticCode(item.code),
                    'reading_frame': ReadingFrame(item.frame),
                    'codon_names': tuple(item.names.values()),
                    }
                for item in self.states.edit.view.iterate()
                if item.split}
            self.update(1, 1, 'text')
            # raise Exception('this crashes the state machine...')
        return self.states.edit.marked.value

    def onFail(self, exception, trace):
        # raise exception
        self.states.wait.logio.writeline('')
        self.states.wait.logio.writeline(trace)
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Critical)
        msgBox.setText(type(exception).__name__)
        msgBox.setInformativeText(str(exception))
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.machine().parent().msgShow(msgBox)

    def skipWait(self):
        # Bypass StepWait since autodetect is disabled
        self.work()
        return True
        skip = self.states.edit.marked.value == 0
        return bool(skip)

    def filterCancel(self, event):
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Question)
        msgBox.setText('Cancel subsetting?')
        msgBox.setStandardButtons(
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.No)
        res = self.machine().parent().msgShow(msgBox)
        return res == QtWidgets.QMessageBox.Yes
