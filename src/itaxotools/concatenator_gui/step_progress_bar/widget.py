# -----------------------------------------------------------------------------
# StepProgressBar - A simple step progress widget for PySide6
# Copyright (C) 2021-2023  Patmanidis Stefanos
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

"""StepProgressBar widget definition"""


from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui

from dataclasses import dataclass, field

from . import states
from . import palette


@palette.palette
class Palette():
    """Theme-aware three-color palette"""

    def bold(self):
        """For active parts"""
        qt_palette = QtGui.QGuiApplication.palette()
        color = qt_palette.color(QtGui.QPalette.Shadow)
        return color

    def base(self):
        """For inactive parts"""
        qt_palette = QtGui.QGuiApplication.palette()
        color = qt_palette.color(QtGui.QPalette.Dark)
        return color

    def weak(self):
        """For background parts"""
        qt_palette = QtGui.QGuiApplication.palette()
        color = qt_palette.color(QtGui.QPalette.Dark).lighter(140)
        return color


@dataclass
class Step():
    """
    Holds information for each step:
    - `text` will be visible on the bar
    - `weight` affects the length of the line after the step
    - `visible` determines if the state will be displayed
    The following are used internally by StepProgressBar:
    - `status` affects text and indicator style
    - `width` is the minimum required for the text
    - `pos` holds the step horizontal position
    """
    text: str = ''
    weight: int = 1
    visible: bool = True
    status: states.AbstractStatus = states.Pending
    width: int = field(repr=False, default=0)
    pos: int = field(repr=False, default=0)

    def drawText(self, painter, palette):
        self.status.drawText(painter, palette, self.text)

    def drawIndicator(self, painter, palette):
        self.status.drawIndicator(painter, palette)


class StepProgressBar(QtWidgets.QWidget):
    """
    Shows user progress in a series of steps towards a goal.
    The active step is highlighted as either active, ongoing or failed.

    Text font, color style and spacing are configurable.
    Ongoing status is animated.

    Usage example
    -------------
    import stepprogressbar
    stepProgressBar = stepprogressbar.StepProgressBar()
    stepProgressBar.addStep('a', 'A', 1)
    stepProgressBar.addStep('b', 'B', 2, False)
    stepProgressBar.addStep('c', 'C')
    stepProgressBar.activateKey('c')
    stepProgressBar.setOngoing()
    """

    def __init__(self, steps=[], font=None, *args, **kwargs):
        """Steps and font may be set on construction"""
        super().__init__(*args, **kwargs)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Minimum,
            QtWidgets.QSizePolicy.Policy.Minimum)
        self.palette = Palette()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.handleTimer)
        self.timerStep = 10
        self.textPadding = 20
        self.verticalPadding = 4
        self.indicatorPadding = 6
        self.steps = []
        self.keys = {}
        self.active = -1
        self.font = font if font is not None else QtGui.QGuiApplication.font()
        self.setSteps(steps)

    @property
    def font(self):
        return self._font

    @font.setter
    def font(self, font):
        self._font = font
        self.metrics = QtGui.QFontMetrics(self.font)
        for step in self.steps:
            self.updateStepWidth(step)

    def updateStepWidth(self, step):
        step.width = self.metrics.horizontalAdvance(step.text)

    def addStep(self, key=None, text='', weight=1, visible=True, step=None):
        """Steps may be added along with a key for faster reference."""
        if step is None:
            step = Step(text=text, weight=weight, visible=visible)
        self.updateStepWidth(step)
        self.steps.append(step)
        if key is not None:
            self.keys[key] = len(self.steps) - 1

    def setSteps(self, strings):
        self.steps = []
        self.keys = {}
        for string in strings:
            self.addStep(string)

    def getVisibleActiveStep(self):
        prevs = [step for step in self.steps[:self.active+1] if step.visible]
        return prevs[-1] if prevs else None

    def activateKey(self, key):
        """Activates the Step identified by given key"""
        index = self.keys[key] if key is not None else -1
        self.activateIndex(index)

    def activateIndex(self, index):
        """Activates the Step identified by given index"""
        index = min(max(index, -1), len(self.steps))
        for i in range(index + 1, len(self.steps)):
            self.steps[i].status = states.Pending
        visibleSteps = [step for step in self.steps if step.visible]
        visibleSteps[0].status = states.Milestone
        visibleSteps[-1].status = states.Milestone
        for i in range(0, index):
            self.steps[i].status = states.Complete
        self.active = index
        self.setStatus(states.Active)
        self.repaint()

    def activateNext(self):
        self.activateIndex(self.active + 1)

    def activatePrevious(self):
        self.activateIndex(self.active - 1)

    def activateFirst(self):
        self.activateIndex(-1)

    def activateFinal(self):
        self.activateIndex(len(self.steps))

    def setStatus(self, status=states.Active):
        """Mark current step with given status"""
        index = self.active
        if index == len(self.steps) - 1 and not self.steps[index].visible:
            return
        if index >= 0 and index < len(self.steps):
            if step := self.getVisibleActiveStep():
                step.status = status
                if status == states.Ongoing:
                    self.timer.start(self.timerStep)
                else:
                    self.timer.stop()
        self.repaint()

    def handleTimer(self):
        self.repaint()

    def paintEvent(self, event):
        painter = QtGui.QPainter()
        painter.begin(self)
        self.draw(painter)
        painter.end()

    def draw(self, painter):
        visibleSteps = [step for step in self.steps if step.visible]
        width = self.size().width()
        height = self.size().height()

        painter.setRenderHint(QtGui.QPainter.Antialiasing)

        totalStepWeight = sum(step.weight for step in visibleSteps[:-1])
        extraWidth = width - self.minimumWidth()
        extraHeight = height - self.minimumHeight()

        cursor = 0
        for step in visibleSteps:
            cursor += self.textPadding
            cursor += step.width / 2
            step.pos = int(cursor)
            cursor += step.width / 2
            cursor += (step.weight / totalStepWeight) * extraWidth

        textY = extraHeight / 2
        textY += self.verticalPadding
        textY += self.metrics.ascent()
        textY = int(textY)
        lineY = height - extraHeight / 2
        lineY -= self.verticalPadding
        lineY -= states.AbstractStatus.indicatorRadius
        lineY = int(lineY)

        painter.setBrush(QtCore.Qt.NoBrush)
        painter.setPen(QtCore.Qt.NoPen)
        painter.setFont(self.font)

        self.drawStepTexts(visibleSteps, painter, textY)
        self.drawStepIndicators(visibleSteps, painter, lineY)
        self.drawStepLines(visibleSteps, painter, lineY)

    def drawStepTexts(self, steps, painter, textY):
        for step in steps:
            painter.save()
            point = QtCore.QPoint(step.pos, textY)
            painter.translate(point)
            step.drawText(painter, self.palette)
            painter.restore()

    def drawStepIndicators(self, steps, painter, lineY):
        for step in steps:
            painter.save()
            point = QtCore.QPoint(step.pos, lineY)
            painter.translate(point)
            step.drawIndicator(painter, self.palette)
            painter.restore()

    def drawStepLines(self, steps, painter, lineY):
        if self.active >= len(self.steps):
            x = self.size().width()
        elif step := self.getVisibleActiveStep():
            x = step.pos
        else:
            x = 0
        for step1, step2 in zip(steps[:-1], steps[1:]):
            x1 = step1.pos
            x1 += step1.status.indicatorRadius
            x1 += self.indicatorPadding
            x2 = step2.pos
            x2 -= step2.status.indicatorRadius
            x2 -= self.indicatorPadding
            self.drawStepLine(painter, x2 < x, x1, x2, lineY)

    def drawStepLine(self, painter, isComplete, x1, x2, y):
        if isComplete:
            color = self.palette.bold
            pen = QtGui.QPen(color, 2, QtCore.Qt.SolidLine)
            painter.setPen(pen)
        else:
            color = self.palette.base
            pen = QtGui.QPen(color, 1, QtCore.Qt.SolidLine)
            painter.setPen(pen)
        painter.drawLine(x1, y, x2, y)

    def minimumWidth(self):
        visibleSteps = [step for step in self.steps if step.visible]
        width = sum(step.width for step in visibleSteps)
        width += self.textPadding * (len(visibleSteps) + 1)
        return width

    def minimumHeight(self):
        height = self.metrics.height()
        height += states.AbstractStatus.indicatorRadius * 2
        height += 3 * self.verticalPadding
        return height

    def sizeHint(self):
        return QtCore.QSize(self.minimumWidth(), self.minimumHeight())
