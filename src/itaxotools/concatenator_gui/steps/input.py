# -----------------------------------------------------------------------------
# ConcatenatorQt - GUI for Concatenator
# Copyright (C) 2021-2024  Patmanidis Stefanos
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

"""StepInput"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtStateMachine
from PySide6 import QtGui

from pathlib import Path

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from itaxotools.common.threading import WorkerThread

from .. import model
from .. import widgets
from .. import step_progress_bar as spb
from .. import step_state_machine as ssm

from ..file_info import file_info_from_path


class InputFrame(common.widgets.Frame):
    def __init__(self, state, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.state = state
        self.setAcceptDrops(True)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.setDropAction(QtCore.Qt.DropAction.LinkAction)
            event.accept()

    def dropEvent(self, event):
        urls = event.mimeData().urls()
        files = [url.toLocalFile() for url in urls]
        self.state.handleAdd(files=files)


class FileItem(widgets.ModelItem):
    fields = [
        'name',
        'format',
        'samples_len',
        'nucleotides',
        'missing',
        'uniform',
        ]

    def __init__(self, parent, file: model.File):
        super().__init__(parent, file)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled)
        self.samples_len = len(self.model.samples)
        font = self.font(0)
        font.setBold(True)
        self.setFont(0, font)


class CharsetItem(widgets.ModelItem):
    fields = [
        'name',
        'format',
        'samples_len',
        'nucleotides',
        'missing',
        'uniform',
        ]
    format = ''

    def __init__(self, parent, charset: model.Charset):
        super().__init__(parent, charset)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemNeverHasChildren)
        self.samples_len = len(self.model.samples)


class StepInputIdle(QtStateMachine.QState):
    def onEntry(self, event):
        parent = self.parent()
        parent.frame.setEnabled(True)
        parent.overlay.setVisible(False)
        action = ssm.NavigateEvent(event).action()
        backwards = bool(action == ssm.NavigateAction.Back)
        parent.footer.setMode(parent.footer.Mode.Middle, backwards=backwards)
        parent.footer.next.setEnabled(bool(parent.data.files))
        parent.stepProgressBar.setStatus(spb.states.Active)
        parent.refresh_contents()


class StepInputBusy(QtStateMachine.QState):
    def onEntry(self, event):
        parent = self.parent()
        parent.frame.setEnabled(False)
        parent.overlay.setVisible(True)
        parent.footer.setMode(parent.footer.Mode.Wait)
        parent.stepProgressBar.setStatus(spb.states.Ongoing)
        parent.spin.start()
        parent.worker.start()

    def onExit(self, event):
        parent = self.parent()
        parent.spin.stop()
        parent.worker.terminate()


class StepInput(ssm.StepState):

    title = 'Import Input Files'
    description = 'Import and inspect sequence files'

    signalAdd = QtCore.Signal()
    signalDone = QtCore.Signal()
    signalRefresh = QtCore.Signal()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.worker = WorkerThread(self.work)
        self.worker.done.connect(self.onDone)
        self.worker.fail.connect(self.onFail)
        self.worker.cancel.connect(self.onCancel)
        self.machine().installWorker(self.worker)
        self.data = model.Concatenation()
        self.files_queue = []
        self.files_ready = []

        self.clear()

    def onDone(self, result):
        while self.files_ready:
            file = self.files_ready.pop()
            if file.path in self.data.files:
                item = self.view.findItems(
                    file.name, QtCore.Qt.MatchExactly, 0)[0]
                self.remove_item(item)
            self.data.files[file.path] = file
            item = FileItem(self.view, file)
            for charset in file.charsets.values():
                CharsetItem(item, charset)
        self.timestamp_set()
        self.signalDone.emit()

    def onFail(self, exception):
        # raise exception
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Critical)
        msgBox.setText('Import failed:')
        msgBox.setInformativeText(str(exception))
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.machine().parent().msgShow(msgBox)
        self.files_queue = []
        self.files_ready = []
        self.signalDone.emit()

    def onCancel(self, exception):
        self.files_queue = []
        self.files_ready = []
        self.signalDone.emit()

    def work(self):
        paths = (f for f in self.data.files if f not in self.files_queue)
        charsets = {
            k for path in paths
            for k in self.data.files[path].charsets.keys()}
        while self.files_queue:
            path = self.files_queue.pop()
            file = file_info_from_path(
                path, self.data.samples, self.worker.check)
            if any(x in charsets for x in file.charsets):
                raise NotImplementedError(f'Duplicate charsets from "{path}"')
            charsets.update(file.charsets.keys())
            self.files_ready.append(file)
            self.worker.check()

    def clear(self):
        self.data = model.Concatenation()
        self.files_queue = []
        self.view.clear()
        self.files.setValue()
        self.charsets.setValue()
        self.samples.setValue()
        self.nucleotides.setValue()

    def draw(self):
        widget = QtWidgets.QWidget()

        text = (
            'Add sequence files by drag-and-drop or by clicking "Import". '
            'Double-click them to inspect their contents.')

        label = QtWidgets.QLabel(text)
        label.setWordWrap(True)

        frame = self.draw_frame()
        overlay = self.draw_overlay()

        stack = QtWidgets.QStackedLayout()
        stack.setStackingMode(QtWidgets.QStackedLayout.StackAll)
        stack.addWidget(frame)
        stack.addWidget(overlay)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(label)
        layout.addLayout(stack, 1)
        layout.setSpacing(24)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)

        return widget

    def draw_summary(self):
        files = widgets.InfoLabel('Files')
        genes = widgets.InfoLabel('Markers')
        samples = widgets.InfoLabel('Samples')
        nucleotides = widgets.InfoLabel('Nucleotides')

        files.setToolTip(
            'The total number of imported files.')
        genes.setToolTip(
            'The total number of markers contained in all files.')
        samples.setToolTip(
            'The number of unique samples across all genes.')
        nucleotides.setToolTip(
            'The total number of nucleotide characters across all files.')

        summary = QtWidgets.QHBoxLayout()
        summary.addWidget(files)
        summary.addWidget(genes)
        summary.addWidget(samples)
        summary.addWidget(nucleotides)
        summary.addStretch(1)
        summary.setSpacing(24)
        summary.setContentsMargins(4, 0, 4, 0)

        self.files = files
        self.charsets = genes
        self.samples = samples
        self.nucleotides = nucleotides

        return summary

    def draw_frame(self):
        view = widgets.TreeWidget()
        view.itemSelectionChanged.connect(self.handleItemSelectionChanged)
        view.setIndentation(12)
        view.setColumnCount(6, 2)
        view.setHeaderLabels([
            'Name', 'Format', 'Samples',
            'Nucleotides', 'Missing', 'Uniform'])

        headerItem = view.headerItem()
        headerItem.setToolTip(0, 'Input filename or gene name')
        headerItem.setToolTip(1, 'Sequence file format')
        headerItem.setToolTip(2, 'Total number of sequences')
        headerItem.setToolTip(3, 'Total number of nucleotide characters')
        headerItem.setToolTip(4, 'Proportion of missing data')
        headerItem.setToolTip(5, 'Are all sequences of the same length?')

        frame = InputFrame(self)

        add = common.widgets.PushButton('&Import', onclick=self.handleAdd)
        remove = common.widgets.PushButton('&Remove')
        remove.clicked.connect(self.handleRemove)
        remove.setEnabled(False)

        QtGui.QShortcut(QtGui.QKeySequence.Delete, view, self.handleRemove)

        search = widgets.ViewSearchWidget(self, view)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(add)
        controls.addWidget(remove)
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
        self.remove = remove
        self.search = search

        return frame

    def draw_overlay(self):
        overlay = common.widgets.Frame()
        overlay.setStyleSheet("""
            Frame {
                background: qlineargradient(x1: 0, y1: -20, x2: 0, y2: 10,
                    stop: 0 transparent, stop: 1 Palette(Window));
                border: 1px solid Palette(Midlight);
            }""")

        banner = common.widgets.Frame()
        banner.setStyleSheet("""
            Frame {
                background: qlineargradient(x1: 0, y1: -20, x2: 0, y2: 10,
                    stop: 0 Palette(Light), stop: 1 Palette(Midlight));
                border: 1px solid Palette(Mid);
            }""")
        banner.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Maximum,
            QtWidgets.QSizePolicy.Policy.Maximum)

        wait = QtWidgets.QLabel('Loading files, please wait...')
        wait.setStyleSheet("""
            color: palette(Text);
            font-size: 12px;
            letter-spacing: 1px;
            """)
        spin = widgets.SpinningCircle()

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(wait)
        layout.addWidget(spin)
        layout.setSpacing(16)
        layout.setContentsMargins(48, 24, 48, 24)
        banner.setLayout(layout)

        layout = QtWidgets.QGridLayout()
        layout.addWidget(banner, 1, 1)
        layout.setColumnStretch(0, 1)
        layout.setColumnStretch(2, 1)
        layout.setRowStretch(0, 1)
        layout.setRowStretch(2, 1)
        overlay.setLayout(layout)

        self.overlay = overlay
        self.spin = spin

        return overlay

    def cog(self):
        super().cog()

        self.states = ssm.AttrDict()
        self.states.idle = StepInputIdle(self)
        self.states.busy = StepInputBusy(self)
        self.setInitialState(self.states.idle)

        transition = QtStateMachine.QSignalTransition(self.signalAdd)
        transition.setTargetState(self.states.busy)
        self.states.idle.addTransition(transition)
        self.transitions.add = transition

        transition = QtStateMachine.QSignalTransition(self.signalRefresh)
        transition.setTargetState(self.states.idle)
        self.states.idle.addTransition(transition)
        self.transitions.refresh = transition

        transition = QtStateMachine.QSignalTransition(self.signalDone)
        transition.setTargetState(self.states.idle)
        self.states.busy.addTransition(transition)
        self.transitions.done = transition

        transition = self.machine().navigateTransition(
            ssm.NavigateAction.Cancel)
        transition.setTargetState(self.states.idle)
        self.states.busy.addTransition(transition)
        self.transitions.cancel = transition

    def refresh_contents(self):
        files = self.data.files.values()
        nucleotides = sum([file.nucleotides for file in files])
        samples = model.DataGroup(self.data.samples)
        samples.merge([file.samples for file in files])
        self.files.setValue(len(self.data.files))
        self.nucleotides.setValue(nucleotides)
        self.charsets.setValue(len(self.data.charsets))
        self.samples.setValue(len(samples))
        self.view.resizeColumnsToContents()

    def remove_item(self, item):
        self.data.remove_file(item.model)
        index = self.view.indexOfTopLevelItem(item)
        if not index < 0:
            self.view.takeTopLevelItem(index)

    def handleAdd(self, checked=False, files=[]):
        if not files:
            (files, _) = QtWidgets.QFileDialog.getOpenFileNames(
                self.machine().parent(),
                self.machine().parent().title + ' - Open File',
                QtCore.QDir.currentPath(),
                'All Files (*)')
        if not files:
            return
        paths = [Path(file) for file in files]
        self.files_queue.extend(paths)
        self.signalAdd.emit()

    def handleRemove(self, checked=False):
        items = [item for item in self.view.selectedItems()
                 if isinstance(item, FileItem)]
        for item in items:
            self.remove_item(item)
        self.timestamp_set()
        self.signalRefresh.emit()

    def handleItemSelectionChanged(self):
        items = self.view.selectedItems()
        enabled = False
        for item in items:
            if isinstance(item, FileItem):
                enabled = True
                break
        self.remove.setEnabled(enabled)
