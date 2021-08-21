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

"""Status indicators for StepProgressBar Steps"""

from PySide6 import QtCore


class AbstractIndicator():

    radius = 10

    @classmethod
    def draw(cls, painter, point):
        raise NotImplementedError()


class Active(AbstractIndicator):

    radius = 6
    radiusInner = 3

    @classmethod
    def draw(cls, painter, point):
        painter.save()
        painter.setPen(QtCore.Qt.NoPen)
        painter.drawEllipse(point, cls.radiusInner, cls.radiusInner)
        painter.restore()
        painter.save()
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.drawEllipse(point, cls.radius, cls.radius)
        painter.restore()


class Pending(AbstractIndicator):
    radius = 4

    @classmethod
    def draw(cls, painter, point):
        painter.save()
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.drawEllipse(point, cls.radius, cls.radius)
        painter.restore()


class Complete(AbstractIndicator):
    pass


class Failed(AbstractIndicator):
    pass


class Ongoing(AbstractIndicator):
    pass
