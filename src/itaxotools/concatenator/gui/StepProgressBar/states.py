# -----------------------------------------------------------------------------
# StepProgressBar - A custom widget for PySide6
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

"""Status classes for StepProgressBar Steps"""

from PySide6 import QtCore


class AbstractStatus():

    indicatorRadius = 10

    @classmethod
    def drawIndicator(cls, painter, point):
        raise NotImplementedError()


class Active(AbstractStatus):

    indicatorRadius = 6
    indicatorInner = 3

    @classmethod
    def drawText(cls, painter, point, text):
        painter.save()
        font = painter.font()
        size = font.pointSize()
        font.setBold(True)
        font.setUnderline(True)
        font.setPointSize(size + 1)
        painter.setFont(font)
        metrics = painter.fontMetrics()
        width = metrics.horizontalAdvance(text)
        descent = metrics.descent()
        painter.drawText(point.x() - width/2, point.y() - descent, text)
        painter.restore()

    @classmethod
    def drawIndicator(cls, painter, point):
        painter.save()
        painter.setPen(QtCore.Qt.NoPen)
        painter.drawEllipse(point, cls.indicatorInner, cls.indicatorInner)
        painter.restore()
        painter.save()
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.drawEllipse(point, cls.indicatorRadius, cls.indicatorRadius)
        painter.restore()


class Pending(AbstractStatus):
    indicatorRadius = 4

    @classmethod
    def drawText(cls, painter, point, text):
        painter.save()
        metrics = painter.fontMetrics()
        width = metrics.horizontalAdvance(text)
        descent = metrics.descent()
        painter.drawText(point.x() - width/2, point.y() - descent, text)
        painter.restore()

    @classmethod
    def drawIndicator(cls, painter, point):
        painter.save()
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.drawEllipse(point, cls.indicatorRadius, cls.indicatorRadius)
        painter.restore()


class Complete(AbstractStatus):
    pass


class Failed(AbstractStatus):
    pass


class Ongoing(AbstractStatus):
    pass
