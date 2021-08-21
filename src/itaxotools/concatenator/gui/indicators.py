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
from PySide6 import QtWidgets
from PySide6 import QtGui


class AbstractIndicator():

    radius = 20

    @classmethod
    def draw(cls):
        raise NotImplementedError()

class Current(AbstractIndicator):
    pass

class Pending(AbstractIndicator):
    pass

class Complete(AbstractIndicator):
    pass

class Failed(AbstractIndicator):
    pass

class Ongoing(AbstractIndicator):
    pass
