# -----------------------------------------------------------------------------
# StepProgressBar - A simple step progress widget for PySide6
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
    """Abstract class for Step status"""

    # This abstract radius is used as a maximum value when drawing
    indicatorRadius = 10

    @classmethod
    def drawText(cls, painter, palette, text):
        """
        Stylized drawing of text for given status.
        The painter must be already centered to the drawing point.
        Painter must have NoPen, NoBrush and a set font.
        May change the painter properties.
        """
        raise NotImplementedError()

    @classmethod
    def drawIndicator(cls, painter, palette):
        """
        Stylized drawing of the graphics indicator for given status.
        The painter must be already centered to the drawing point,
        with the vertical position at the desired text baseline.
        Painter must have NoPen and NoBrush.
        May change the painter properties.
        """
        raise NotImplementedError()


class Pending(AbstractStatus):

    indicatorRadius = 3

    @classmethod
    def drawText(cls, painter, palette, text):
        painter.setPen(QtGui.QPen(palette.base, 0, QtCore.Qt.SolidLine))
        width = painter.fontMetrics().horizontalAdvance(text)
        painter.drawText(-width/2, 0, text)

    @classmethod
    def drawIndicator(cls, painter, palette):
        point = QtCore.QPoint(0, 0)
        painter.setPen(QtGui.QPen(palette.weak, 2, QtCore.Qt.SolidLine))
        painter.drawEllipse(point, cls.indicatorRadius, cls.indicatorRadius)


class Final(Pending):

    indicatorRadius = 4


class Complete(Pending):

    indicatorRadius = 6

    @classmethod
    def drawIndicator(cls, painter, palette):
        """Draw a check mark"""
        polygon = QtGui.QPolygon([
            QtCore.QPoint(-4, 0),
            QtCore.QPoint(-2, 3),
            QtCore.QPoint(6, -6)])
        painter.setPen(QtGui.QPen(palette.bold, 2, QtCore.Qt.SolidLine))
        painter.drawPolyline(polygon)


class Active(AbstractStatus):

    indicatorRadius = 6
    indicatorInner = 3

    @classmethod
    def drawText(cls, painter, palette, text):
        """Highlight text"""
        painter.setPen(QtGui.QPen(palette.bold, 2, QtCore.Qt.SolidLine))
        font = painter.font()
        size = font.pointSize()
        font.setBold(True)
        font.setUnderline(True)
        # font.setPointSize(size + 1)
        painter.setFont(font)
        width = painter.fontMetrics().horizontalAdvance(text)
        painter.drawText(-width/2, 0, text)

    @classmethod
    def drawIndicator(cls, painter, palette):
        """Large dotted circle"""
        point = QtCore.QPoint(0, 0)

        painter.setPen(QtCore.Qt.NoPen)
        painter.setBrush(palette.bold)
        painter.drawEllipse(point, cls.indicatorInner, cls.indicatorInner)

        painter.setPen(QtGui.QPen(palette.bold, 2, QtCore.Qt.SolidLine))
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.drawEllipse(point, cls.indicatorRadius, cls.indicatorRadius)


class Failed(Active):

    indicatorRadius = 6
    indicatorInner = 4

    @classmethod
    def drawIndicator(cls, painter, palette):
        """X-mark"""
        rad = cls.indicatorInner
        painter.setPen(QtGui.QPen(palette.bold, 3, QtCore.Qt.SolidLine))
        painter.drawLine(-rad, -rad, rad, rad)
        painter.drawLine(-rad, rad, rad, -rad)


class Ongoing(Active):

    indicatorRadius = 6
    indicatorSpan = 120
    indicatorPeriod = 2

    @classmethod
    def drawIndicator(cls, painter, palette):
        """Spinning circle"""
        rad = cls.indicatorRadius
        rect = QtCore.QRect(-rad, -rad, 2 * rad, 2 * rad)

        painter.setPen(QtGui.QPen(palette.weak, 2, QtCore.Qt.SolidLine))
        painter.drawEllipse(rect)

        period_ns = int(cls.indicatorPeriod * 10**9)
        ns = time_ns() % period_ns
        degrees = - 360 * ns / period_ns
        painter.setPen(QtGui.QPen(palette.bold, 2, QtCore.Qt.SolidLine))
        painter.drawArc(rect, degrees * 16, cls.indicatorSpan * 16)
