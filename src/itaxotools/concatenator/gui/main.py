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

"""Long module description"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtStateMachine
from PySide6 import QtGui

import tempfile
import pathlib
import shutil
import enum
import sys

from lorem_text import lorem
from random import randint
from time import time_ns

from itaxotools.common import utility
from itaxotools.common import widgets
from itaxotools.common import resources
from itaxotools.common import io

from . import step_progress_bar as spb
from . import step_state_machine as ssm


try:
    import importlib.resources
    _resource_path = importlib.resources.files(resources)
    _package_path = importlib.resources.files(__package__)
except Exception:
    if hasattr(sys, '_MEIPASS'):
        _resource_path = (pathlib.Path(sys._MEIPASS) /
            'itaxotools' / 'common' / 'resources')
        _package_path = (pathlib.Path(sys._MEIPASS) /
            'itaxotools' / 'fasttreepy' / 'gui')
    else:
        import os
        _resource_path = pathlib.Path(os.path.dirname(resources.__file__))
        _package_path = pathlib.Path(os.path.dirname(__file__))

def get_resource(path):
    return str(_resource_path / path)
def get_icon(path):
    return str(_resource_path / 'icons/svg' / path)


def dummy_work(self, count, max, lines, period):
    with io.redirect(sys, 'stdout', self.logio):
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


class SpinningCircle(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.handleTimer)
        self.timerStep = 10
        self.radius = 8
        self.period = 2
        self.span = 120
        self.width = 2

    def start(self):
        self.timer.start(self.timerStep)

    def stop(self):
        self.timer.stop()

    def handleTimer(self):
        self.repaint()

    def sizeHint(self):
        diameter = (self.radius + self.width) * 2
        return QtCore.QSize(diameter, diameter)

    def paintEvent(self, event):
        painter = QtGui.QPainter()
        painter.begin(self)

        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setBrush(QtCore.Qt.NoBrush)

        x = self.size().width()/2
        y = self.size().height()/2
        painter.translate(QtCore.QPoint(x, y))

        palette = QtGui.QGuiApplication.palette()
        weak = palette.color(QtGui.QPalette.Mid)
        bold = palette.color(QtGui.QPalette.Shadow)

        rad = self.radius
        rect = QtCore.QRect(-rad, -rad, 2 * rad, 2 * rad)

        painter.setPen(QtGui.QPen(weak, self.width, QtCore.Qt.SolidLine))
        painter.drawEllipse(rect)

        period_ns = int(self.period * 10**9)
        ns = time_ns() % period_ns
        degrees = - 360 * ns / period_ns
        painter.setPen(QtGui.QPen(bold, self.width, QtCore.Qt.SolidLine))
        painter.drawArc(rect, degrees * 16, self.span * 16)

        painter.end()


class InfoLabel(QtWidgets.QLabel):
    def __init__(self, text):
        super().__init__()
        self.prefix = text
        self.setValue('-')

    def setValue(self, value):
        self.value = value
        if isinstance(value, int):
            value = f'{value:,}'
        self.setText(f'{self.prefix}: {value}')

class StepAbout(ssm.StepState):
    pass


class StepInput(ssm.StepState):

    title = 'Select Input Files'
    description = 'Add and inspect sequence files'

    signalAdd = QtCore.Signal()
    signalDone = QtCore.Signal()
    signalRefresh = QtCore.Signal()

    class StepInputIdle(QtStateMachine.QState):
        def onEntry(self, event):
            parent = self.parent()
            parent.frame.setEnabled(True)
            parent.overlay.setVisible(False)
            parent.footer.setMode(parent.footer.Mode.Middle)
            parent.footer.next.setEnabled(bool(parent.data.files))
            parent.stepProgressBar.setStatus(spb.states.Active)

    class StepInputBusy(QtStateMachine.QState):
        def onEntry(self, event):
            parent = self.parent()
            parent.frame.setEnabled(False)
            parent.overlay.setVisible(True)
            parent.footer.setMode(parent.footer.Mode.Wait)
            parent.stepProgressBar.setStatus(spb.states.Ongoing)

            parent.worker.start()

        def onExit(self, event):
            parent = self.parent()
            parent.worker.terminate()

    class InputFrame(widgets.Frame):
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

    class FileItem(QtWidgets.QTreeWidgetItem):
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
            self.file = None
            for col in range(1, 6):
                self.setTextAlignment(col, QtCore.Qt.AlignRight)
            self.name = lorem.words(randint(2, 6)).replace(' ', '_')
            self.format = ['Fasta', 'Phylip', 'Nexus'][randint(0,2)]
            self.samples = randint(6,300)
            self.nucleotides = randint(200, 3000000)
            self.uniform = ['Yes', 'No'][randint(0,1)]
            self.missing = f'{randint(0,99)}%'

        def __setattr__(self, attr, value):
            super().__setattr__(attr, value)
            if attr in self.map:
                if isinstance(value, int):
                    value = f'{value:,}'
                self.setText(self.map[attr], str(value))

    class SetItem(QtWidgets.QTreeWidgetItem):
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
            for col in range(1, 6):
                self.setTextAlignment(col, QtCore.Qt.AlignRight)
            self.name = lorem.words(randint(2, 6)).replace(' ', '_')
            self.format = '-'
            self.samples = randint(6,300)
            self.nucleotides = randint(200, 3000000)
            self.uniform = ['Yes', 'No'][randint(0,1)]
            self.missing = f'{randint(0,99)}%'

        def __setattr__(self, attr, value):
            super().__setattr__(attr, value)
            if attr in self.map:
                if isinstance(value, int):
                    value = f'{value:,}'
                self.setText(self.map[attr], str(value))

    class DataObject(object):
        pass

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = self.DataObject()
        self.data.files = []
        self.data.files_pending = []
        self.machine().installWorker(self)

    def onDone(self, result):
        total_nucleotides = 0
        total_samples = 0
        total_sets = 0
        while self.data.files_pending:
            file = self.data.files_pending.pop()
            self.data.files.append(file)
            item = self.FileItem(self.view)
            item.file = file
            item.name = file.name
            for i in range(1, randint(1, 99)):
                self.SetItem(item)
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

        self.files.setValue(len(self.data.files))
        self.nucleotides.setValue(total_nucleotides)
        self.samples.setValue(total_samples)
        self.sets.setValue(total_sets)
        self.signalDone.emit()

    def onFail(self, exception):
        self.signalDone.emit()

    def onCancel(self, exception):
        self.signalDone.emit()

    def work(self):
        time = randint(500, 3000)
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

        self.frame.setEnabled(False)
        self.spin.start()

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
        files = InfoLabel('Input Files')
        sets = InfoLabel('Character Sets')
        samples = InfoLabel('Unique Samples')
        nucleotides = InfoLabel('Nucleotides')

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
        frame = self.InputFrame(self)

        add = widgets.PushButton('&Add Files')
        add.clicked.connect(self.handleAdd)
        remove = widgets.PushButton('&Remove Files')
        remove.clicked.connect(self.handleRemove)

        action = QtGui.QAction('Search', self)
        pixmap = widgets.VectorPixmap(get_icon('search.svg'),
            colormap=self.machine().parent().colormap_icon_light)
        action.setIcon(pixmap)
        action.setShortcut(QtGui.QKeySequence.FindNext)
        action.setStatusTip('Search')
        action.triggered.connect(self.handleSearch)

        search = widgets.SearchWidget()
        search.setSearchAction(action)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(add)
        controls.addWidget(remove)
        controls.addStretch(1)
        controls.addWidget(search)
        controls.setSpacing(8)
        controls.setContentsMargins(0, 0, 0, 0)

        view = QtWidgets.QTreeWidget()
        view.itemSelectionChanged.connect(self.handleItemSelectionChanged)
        view.setSelectionMode(QtWidgets.QTreeWidget.ExtendedSelection)
        view.setIndentation(12)
        view.setAlternatingRowColors(True)
        view.setSortingEnabled(True)
        view.setColumnCount(6)
        view.setHeaderLabels([
            ' Name', ' Format', ' Samples',
            ' Nucleotides', ' Uniform', ' Missing'])

        headerItem = view.headerItem()
        for col in range(1, 6):
            headerItem.setTextAlignment(col, QtCore.Qt.AlignRight)

        summary = self.draw_summary()

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(controls)
        layout.addWidget(view, 1)
        layout.addLayout(summary)
        layout.setSpacing(8)
        layout.setContentsMargins(16, 16, 16, 8)
        frame.setLayout(layout)

        self.view = view
        self.frame = frame
        self.remove = remove
        self.search = search

        return frame

    def draw_overlay(self):
        overlay = widgets.Frame()
        overlay.setStyleSheet("""
            Frame {
                background: qlineargradient(x1: 0, y1: -20, x2: 0, y2: 10,
                    stop: 0 transparent, stop: 1 Palette(Window));
                border: 1px solid Palette(Midlight);
            }""")

        banner = widgets.Frame()
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
        spin = SpinningCircle()

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
        self.states['idle'] = self.StepInputIdle(self)
        self.states['busy'] = self.StepInputBusy(self)
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
                      if isinstance(item, self.FileItem)]
        for item in items:
            self.data.files.remove(item.file)
            index = self.view.indexOfTopLevelItem(item)
            self.view.takeTopLevelItem(index)
        self.signalRefresh.emit()

    def handleSearch(self, checked=False):
        print('search')

    def handleItemSelectionChanged(self):
        items = self.view.selectedItems()
        enabled = False
        for item in items:
            if isinstance(item, self.FileItem):
                enabled = True
                break
        self.remove.setEnabled(enabled)

class StepFilter(ssm.StepState):
    pass


class StepAlign(ssm.StepTriState):
    title = 'Sequence Alignment'

    class StepEdit(ssm.StepSubState):
        description = 'Options, options'

    class StepFail(ssm.StepSubState):
        description = 'Task failed'

    class StepWait(ssm.StepSubState):
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
            logger = widgets.TextEditLogger()
            logger.setFont(QtGui.QFontDatabase.systemFont(
                QtGui.QFontDatabase.FixedFont))
            logger.document().setDocumentMargin(10)
            logio = io.TextEditLoggerIO(logger)

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

    def updateProgress(self, tuple):
        x, m, w = tuple
        self.progtext.setText(f'Sequence {x}/{m}: {w}')
        self.progbar.setMaximum(m)
        self.progbar.setValue(x)

    def work(self):
        dummy_work(self, 42, 200, 10, 10)

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


class StepCodons(ssm.StepTriState):
    pass


class StepExport(ssm.StepState):
    pass


class StepDone(ssm.StepState):
    pass


class Main(widgets.ToolDialog):
    """Main window, handles everything"""

    def __init__(self, parent=None, init=None):
        super(Main, self).__init__(parent)

        self.title = 'Concatenator'
        self._temp = None
        self.temp = None

        icon = QtGui.QIcon(get_resource('logos/ico/concatenator.ico'))
        self.setWindowIcon(icon)
        self.setWindowTitle(self.title)
        self.resize(840, 540)

        self.process = None
        self.machine = None
        self.skin()
        self.draw()
        self.act()
        self.cog()

        if init is not None:
            self.machine.started.connect(init)

    def __getstate__(self):
        return (None,)

    def __setstate__(self, state):
        pass

    def skin(self):
        """Configure widget appearance"""
        color = {
            'white':  '#ffffff',
            'light':  '#eff1ee',
            'beige':  '#e1e0de',
            'gray':   '#abaaa8',
            'iron':   '#8b8d8a',
            'black':  '#454241',
            'red':    '#ee4e5f',
            'pink':   '#eb9597',
            'orange': '#eb6a4a',
            'brown':  '#655c5d',
            'green':  '#00ff00',
            }
        # using green for debugging
        palette = QtGui.QGuiApplication.palette()
        scheme = {
            QtGui.QPalette.Active: {
                QtGui.QPalette.Window: 'light',
                QtGui.QPalette.WindowText: 'black',
                QtGui.QPalette.Base: 'white',
                QtGui.QPalette.AlternateBase: 'light',
                QtGui.QPalette.PlaceholderText: 'gray',
                QtGui.QPalette.Text: 'black',
                QtGui.QPalette.Button: 'light',
                QtGui.QPalette.ButtonText: 'black',
                QtGui.QPalette.Light: 'white',
                QtGui.QPalette.Midlight: 'beige',
                QtGui.QPalette.Mid: 'gray',
                QtGui.QPalette.Dark: 'iron',
                QtGui.QPalette.Shadow: 'brown',
                QtGui.QPalette.Highlight: 'red',
                QtGui.QPalette.HighlightedText: 'white',
                # These work on linux only?
                QtGui.QPalette.ToolTipBase: 'beige',
                QtGui.QPalette.ToolTipText: 'brown',
                # These seem bugged anyway
                QtGui.QPalette.BrightText: 'green',
                QtGui.QPalette.Link: 'red',
                QtGui.QPalette.LinkVisited: 'pink',
                },
            QtGui.QPalette.Disabled: {
                QtGui.QPalette.Window: 'light',
                QtGui.QPalette.WindowText: 'iron',
                QtGui.QPalette.Base: 'white',
                QtGui.QPalette.AlternateBase: 'light',
                QtGui.QPalette.PlaceholderText: 'beige',
                QtGui.QPalette.Text: 'iron',
                QtGui.QPalette.Button: 'light',
                QtGui.QPalette.ButtonText: 'gray',
                QtGui.QPalette.Light: 'white',
                QtGui.QPalette.Midlight: 'beige',
                QtGui.QPalette.Mid: 'gray',
                QtGui.QPalette.Dark: 'iron',
                QtGui.QPalette.Shadow: 'brown',
                QtGui.QPalette.Highlight: 'pink',
                QtGui.QPalette.HighlightedText: 'white',
                # These seem bugged anyway
                QtGui.QPalette.BrightText: 'green',
                QtGui.QPalette.ToolTipBase: 'green',
                QtGui.QPalette.ToolTipText: 'green',
                QtGui.QPalette.Link: 'green',
                QtGui.QPalette.LinkVisited: 'green',
                },
            }
        scheme[QtGui.QPalette.Inactive] = scheme[QtGui.QPalette.Active]
        for group in scheme:
            for role in scheme[group]:
                palette.setColor(
                    group, role, QtGui.QColor(color[scheme[group][role]]))
        QtGui.QGuiApplication.setPalette(palette)

        self.colormap = {
            widgets.VectorIcon.Normal: {
                '#000': color['black'],
                '#f00': color['red'],
                },
            widgets.VectorIcon.Disabled: {
                '#000': color['gray'],
                '#f00': color['orange'],
                },
            }
        self.colormap_icon = {
            '#000': color['black'],
            '#f00': color['red'],
            '#f88': color['pink'],
            }
        self.colormap_icon_light = {
            '#000': color['gray'],
            '#ff0000': color['red'],
            '#ffa500': color['pink'],
            }

    def draw(self):
        """Draw all widgets"""

        self.header = widgets.Header()
        self.header.logoTool = widgets.VectorPixmap(
            get_resource('logos/svg/concatenator.svg'),
            size=QtCore.QSize(48, 48), colormap=self.colormap_icon)
        self.header.logoProject = QtGui.QPixmap(
            get_resource('logos/png/itaxotools-micrologo.png'))

        self.header.setTool(
            title='Concatenator',
            citation='by V. Kharchev and S. Patmanidis',
            description='Sequence alignment file manipulation')

        self.stepProgressBar = spb.StepProgressBar()
        font = QtGui.QGuiApplication.font()
        font.setPointSize(9)
        font.setLetterSpacing(QtGui.QFont.AbsoluteSpacing, 1)
        font.setBold(True)
        self.stepProgressBar.font = font

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.stepProgressBar)
        layout.setContentsMargins(0, 8, 0, 8)
        self.header.widget.setLayout(layout)

        self.body = QtWidgets.QStackedLayout()
        self.footer = widgets.NavigationFooter()

        layout = QtWidgets.QGridLayout()
        layout.addWidget(self.header, 0, 0, 1, 5)
        layout.addLayout(self.body, 2, 2)
        layout.addWidget(self.footer, 4, 1, 1, 3)
        layout.setRowStretch(2, 1)
        layout.setColumnMinimumWidth(0, 8)
        layout.setColumnMinimumWidth(4, 8)
        layout.setColumnMinimumWidth(1, 32)
        layout.setColumnMinimumWidth(3, 32)
        layout.setRowMinimumHeight(1, 24)
        layout.setRowMinimumHeight(3, 24)
        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

    def act(self):
        pass

    def cog(self):
        """Initialize state machine"""

        self.machine = ssm.StepStateMachine(
            self, self.stepProgressBar, self.header, self.footer, self.body)
        self.machine.addStep('about', 'About', 1, False, StepAbout)
        self.machine.addStep('input', 'Input', 1, True, StepInput)
        self.machine.addStep('filter', 'Filter', 1, True, StepFilter)
        self.machine.addStep('align', 'Align', 1, True, StepAlign)
        self.machine.addStep('codons', 'Codons', 1, True, StepCodons)
        self.machine.addStep('export', 'Export', 1, True, StepExport)
        self.machine.addStep('done', 'Done', 1, False, StepDone)

        # self.machine.setInitialState(self.machine.states['align'])

        self.machine.start()

    def onReject(self):
        return self.machine.cancel()

    def postReject(self):
        self.machine.terminate()
