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

"""Diagnostic record system"""


from typing import Iterator, Optional
from enum import Enum, auto
from pathlib import Path

from PySide6.QtCore import QAbstractItemModel
from PySide6.QtWidgets import QWidget


class RecordFlag(Enum):
    Info = auto()
    Warn = auto()
    Fail = auto()


class RecordData:
    def __init__(self, data: object):
        self.data = data

    def export(self, path: Path) -> None:
        raise NotImplementedError()

    def model(self) -> QAbstractItemModel:
        raise NotImplementedError()

    def view(self) -> QWidget:
        raise NotImplementedError()


class Record:
    def __init__(
        self,
        type: RecordFlag,
        title: str,
        description: str = '',
        data: Optional[RecordData] = None,
    ):
        self.type = type
        self.title = title
        self.description = description
        self.data = data

    def __str__(self):
        return f"Record({self.type.name}: {repr(self.title)})"


class RecordLog:
    def __init__(self, records: Optional[Iterator[Record]] = None):
        self.records = list()
        if records is not None:
            self.records = list(records)

    def __iter__(self):
        return iter(self.records)

    def __str__(self):
        return f"RecordLog({', '.join(str(x) for x in self)})"

    def append(self, record: Optional[Record]) -> None:
        if record is not None:
            self.records.append(record)
