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

"""StepInput"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtStateMachine

import pathlib

from lorem_text import lorem
from random import randint

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from .. import widgets
from .. import step_progress_bar as spb
from .. import step_state_machine as ssm


class StepInputIdle(QtStateMachine.QState):
    def onEntry(self, event):
        parent = self.parent()
        parent.frame.setEnabled(True)
        parent.overlay.setVisible(False)
        parent.footer.setMode(parent.footer.Mode.Middle)
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
        filenames = [url.toLocalFile() for url in urls]
        self.state.handleAdd(filenames=filenames)


class FileItem(widgets.WidgetItem):
    map = {
        'name': 0,
        'format': 1,
        'samples': 2,
        'nucleotides': 3,
        'uniform': 4,
        'missing': 5,
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled)
        self.file = None
        alignment = QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter
        for col in range(1, 6):
            self.setTextAlignment(col, alignment)
        font = self.font(0)
        font.setBold(True)
        self.setFont(0, font)
        self.name = lorem.words(randint(2, 6)).replace(' ', '_')
        self.format = ['Fasta', 'Phylip', 'Nexus'][randint(0, 2)]
        self.samples = randint(6, 300)
        self.nucleotides = randint(200, 3000000)
        self.uniform = ['Yes', 'No'][randint(0, 1)]
        self.missing = randint(0, 9999) / 10000


class SetItem(widgets.WidgetItem):
    map = {
        'name': 0,
        'format': 1,
        'samples': 2,
        'nucleotides': 3,
        'uniform': 4,
        'missing': 5,
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemNeverHasChildren)
        alignment = QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter
        for col in range(1, 6):
            self.setTextAlignment(col, alignment)
        self.name = lorem.words(randint(2, 6)).replace(' ', '_')
        self.format = '-'
        self.samples = randint(6, 300)
        self.nucleotides = randint(200, 3000000)
        self.uniform = ['Yes', 'No'][randint(0, 1)]
        self.missing = randint(0, 9999) / 10000


class DataObject(object):
    pass


class StepInput(ssm.StepState):

    title = 'Select Input Files'
    description = 'Add and inspect sequence files'

    signalAdd = QtCore.Signal()
    signalDone = QtCore.Signal()
    signalRefresh = QtCore.Signal()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = DataObject()
        self.data.files = []
        self.data.files_pending = []
        self.machine().installWorker(self)

    def onDone(self, result):
        total_nucleotides = 0
        total_samples = 0
        total_sets = 0
        while self.data.files_pending:
            file = self.data.files_pending.pop()
            if file in self.data.files:
                continue
            self.data.files.append(file)
            item = FileItem(self.view)
            item.file = file
            item.name = file.name
            for i in range(1, randint(100, 990)):
                SetItem(item)
            nucleotides = 0
            samples = item.samples
            for i in range(0, item.childCount()):
                nucleotides += item.child(i).nucleotides
                samples = max(samples, item.child(i).samples)
                total_sets += 1
            item.nucleotides = nucleotides
            item.samples = samples
            total_nucleotides += nucleotides
            total_samples += samples
        self.nucleotides.setValue(total_nucleotides)
        self.samples.setValue(total_samples)
        self.sets.setValue(total_sets)

        self.signalDone.emit()

    def onFail(self, exception):
        self.signalDone.emit()

    def onCancel(self, exception):
        self.signalDone.emit()

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
        overlay = self.draw_overlay()

        # effect = QtWidgets.QGraphicsOpacityEffect(self)
        # effect.setOpacity(0.7)
        # overlay.setGraphicsEffect(effect)

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
        files = widgets.InfoLabel('Input Files')
        sets = widgets.InfoLabel('Character Sets')
        samples = widgets.InfoLabel('Unique Samples')
        nucleotides = widgets.InfoLabel('Nucleotides')

        for item in [files, sets, samples, nucleotides]:
            item.setToolTip(lorem.words(randint(5, 15)))

        summary = QtWidgets.QHBoxLayout()
        summary.addWidget(files)
        summary.addWidget(sets)
        summary.addWidget(samples)
        summary.addWidget(nucleotides)
        summary.addStretch(1)
        summary.setSpacing(24)
        summary.setContentsMargins(4, 0, 4, 0)

        self.files = files
        self.sets = sets
        self.samples = samples
        self.nucleotides = nucleotides

        return summary

    def draw_frame(self):
        view = widgets.TreeWidget()
        view.itemSelectionChanged.connect(self.handleItemSelectionChanged)
        view.setIndentation(12)
        view.setColumnCount(6)
        view.setHeaderLabels([
            ' Name', ' Format', ' Samples',
            ' Nucleotides', ' Uniform', ' Missing'])

        headerItem = view.headerItem()
        headerItem.setToolTip(0, lorem.words(13))
        for col in range(1, 6):
            headerItem.setToolTip(col, lorem.words(randint(5, 15)))

        frame = InputFrame(self)

        add = common.widgets.PushButton('&Add Files')
        add.clicked.connect(self.handleAdd)
        remove = common.widgets.PushButton('&Remove Files')
        remove.clicked.connect(self.handleRemove)
        remove.setEnabled(False)

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

        self.states = dict()
        self.states['idle'] = StepInputIdle(self)
        self.states['busy'] = StepInputBusy(self)
        self.setInitialState(self.states['idle'])

        transition = QtStateMachine.QSignalTransition(self.signalAdd)
        transition.setTargetState(self.states['busy'])
        self.states['idle'].addTransition(transition)
        self.transitions['add'] = transition

        transition = QtStateMachine.QSignalTransition(self.signalRefresh)
        transition.setTargetState(self.states['idle'])
        self.states['idle'].addTransition(transition)
        self.transitions['refresh'] = transition

        transition = QtStateMachine.QSignalTransition(self.signalDone)
        transition.setTargetState(self.states['idle'])
        self.states['busy'].addTransition(transition)
        self.transitions['done'] = transition

        transition = ssm.NavigateTransition(ssm.NavigateEvent.Event.Cancel)
        transition.setTargetState(self.states['idle'])
        self.states['busy'].addTransition(transition)
        self.transitions['cancel'] = transition

    def refresh_contents(self):
        self.files.setValue(len(self.data.files))
        self.view.resizeColumnsToContents()

    def handleAdd(self, checked=False, filenames=[]):
        if not filenames:
            (filenames, _) = QtWidgets.QFileDialog.getOpenFileNames(
                self.machine().parent(),
                self.machine().parent().title + ' - Open File',
                QtCore.QDir.currentPath(),
                'All Files (*)')
        if not filenames:
            return
        paths = [pathlib.Path(filename) for filename in filenames]
        self.data.files_pending.extend(paths)
        self.signalAdd.emit()

    def handleRemove(self, checked=False):
        items = [item for item in self.view.selectedItems()
                 if isinstance(item, FileItem)]
        for item in items:
            self.data.files.remove(item.file)
            index = self.view.indexOfTopLevelItem(item)
            if not index < 0:
                self.view.takeTopLevelItem(index)
        self.signalRefresh.emit()

    def handleItemSelectionChanged(self):
        items = self.view.selectedItems()
        enabled = False
        for item in items:
            if isinstance(item, FileItem):
                enabled = True
                break
        self.remove.setEnabled(enabled)
