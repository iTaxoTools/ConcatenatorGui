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

from time import time_ns


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
        painter.translate(point)
        palette = QtGui.QGuiApplication.palette()
        color = palette.color(QtGui.QPalette.Dark)
        pen = QtGui.QPen(color, 2, QtCore.Qt.SolidLine)
        painter.setPen(pen)
        metrics = painter.fontMetrics()
        width = metrics.horizontalAdvance(text)
        descent = metrics.descent()
        painter.drawText(-width/2, -descent, text)
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

    indicatorRadius = 4


class Complete(Pending):

    indicatorRadius = 6

    @classmethod
    def drawIndicator(cls, painter, point):
        point1 = QtCore.QPoint(-4, 0)
        point2 = QtCore.QPoint(-2, 3)
        point3 = QtCore.QPoint(6, -6)
        polygon = QtGui.QPolygon([point1, point2, point3])
        painter.save()
        painter.translate(point)
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.drawPolyline(polygon)
        painter.restore()


class Active(AbstractStatus):

    indicatorRadius = 6
    indicatorInner = 3

    @classmethod
    def drawText(cls, painter, point, text):
        painter.save()
        painter.translate(point)
        font = painter.font()
        size = font.pointSize()
        font.setBold(True)
        font.setUnderline(True)
        # font.setPointSize(size + 1)
        painter.setFont(font)
        metrics = painter.fontMetrics()
        width = metrics.horizontalAdvance(text)
        descent = metrics.descent()
        painter.drawText(-width/2, -descent, text)
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
         rad = 4
         point1 = QtCore.QPoint(-rad, -rad)
         point2 = QtCore.QPoint(rad, rad)
         point3 = QtCore.QPoint(-rad, rad)
         point4 = QtCore.QPoint(rad, -rad)
         painter.save()
         painter.translate(point)
         pen = painter.pen()
         pen.setWidth(3)
         painter.setPen(pen)
         painter.setBrush(QtCore.Qt.NoBrush)
         painter.drawLine(point1, point2)
         painter.drawLine(point3, point4)
         painter.restore()


class Ongoing(Active):

    indicatorRadius = 6
    indicatorSpan = 120
    indicatorSpeed = 0.5

    @classmethod
    def drawIndicator(cls, painter, point):
        rad = cls.indicatorRadius
        rect = QtCore.QRect(-rad, -rad, 2 * rad, 2 * rad)
        painter.save()
        painter.translate(point)
        color = painter.pen().color().lighter(200)
        pen = QtGui.QPen(color, 2, QtCore.Qt.SolidLine)
        painter.setPen(pen)
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.drawEllipse(rect)
        painter.restore()
        painter.save()
        painter.translate(point)
        painter.setBrush(QtCore.Qt.NoBrush)
        rect = QtCore.QRect(-rad, -rad, 2*rad, 2*rad)
        ns = (time_ns() * cls.indicatorSpeed) % 10**9
        degrees = - 360 * ns / 10**9
        painter.drawArc(rect, degrees * 16, cls.indicatorSpan * 16)
        painter.restore()
