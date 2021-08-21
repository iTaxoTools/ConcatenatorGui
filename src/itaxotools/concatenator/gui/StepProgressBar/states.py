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
from PySide6 import QtGui


class AbstractStatus():

    indicatorRadius = 10

    @classmethod
    def drawText(cls, painter, point, text):
        raise NotImplementedError()

    @classmethod
    def drawIndicator(cls, painter, point):
        raise NotImplementedError()


class Pending(AbstractStatus):

    indicatorRadius = 3

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
        palette = QtGui.QGuiApplication.palette()
        color = palette.color(QtGui.QPalette.Dark)
        pen = QtGui.QPen(color, 2, QtCore.Qt.SolidLine)
        painter.setPen(pen)
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.drawEllipse(point, cls.indicatorRadius, cls.indicatorRadius)
        painter.restore()


class Final(Pending):

    indicatorRadius = 5


class Complete(Pending):

    indicatorRadius = 6

    @classmethod
    def drawIndicator(cls, painter, point):
        point1 = point + QtCore.QPoint(-4, 0)
        point2 = point + QtCore.QPoint(-2, 3)
        point3 = point + QtCore.QPoint(5, -5)
        polygon = QtGui.QPolygon([point1, point2, point3])
        painter.save()
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.drawPolyline(polygon)
        painter.restore()


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
        font.setPointSize(10)
        painter.setFont(font)
        metrics = painter.fontMetrics()
        width = metrics.horizontalAdvance(text)
        descent = metrics.descent()
        painter.drawText(point.x() - width/2, point.y() - descent - 1, text)
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


class Failed(Active):

     indicatorRadius = 6

     @classmethod
     def drawIndicator(cls, painter, point):
         a = 4
         point1 = point + QtCore.QPoint(-a, -a)
         point2 = point + QtCore.QPoint(a, a)
         point3 = point + QtCore.QPoint(-a, a)
         point4 = point + QtCore.QPoint(a, -a)
         painter.save()
         pen = painter.pen()
         pen.setWidth(3)
         painter.setPen(pen)
         painter.setBrush(QtCore.Qt.NoBrush)
         painter.drawLine(point1, point2)
         painter.drawLine(point3, point4)
         painter.restore()


class Ongoing(AbstractStatus):
    pass
