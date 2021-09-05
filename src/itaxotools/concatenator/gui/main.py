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

"""Main dialog window"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui

from lorem_text import lorem
from random import randint

import sys

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa
import itaxotools.common.io # noqa

from . import step_progress_bar as spb
from . import step_state_machine as ssm

from .steps.input import StepInput


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


class StepAbout(ssm.StepState):
    pass


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


class Main(common.widgets.ToolDialog):
    """Main window, handles everything"""

    def __init__(self, parent=None, init=None):
        super(Main, self).__init__(parent)

        self.title = 'Concatenator'
        self._temp = None
        self.temp = None

        icon = QtGui.QIcon(common.resources.get('logos/ico/concatenator.ico'))
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
            common.widgets.VectorIcon.Normal: {
                '#000': color['black'],
                '#f00': color['red'],
                },
            common.widgets.VectorIcon.Disabled: {
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

        self.header = common.widgets.Header()
        self.header.logoTool = common.widgets.VectorPixmap(
            common.resources.get('logos/svg/concatenator.svg'),
            size=QtCore.QSize(48, 48), colormap=self.colormap_icon)
        self.header.logoProject = QtGui.QPixmap(
            common.resources.get('logos/png/itaxotools-logo-64px.png'))

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
        self.footer = common.widgets.NavigationFooter()

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
