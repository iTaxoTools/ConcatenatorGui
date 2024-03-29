# -----------------------------------------------------------------------------
# ConcatenatorQt - GUI for Concatenator
# Copyright (C) 2021-2024  Patmanidis Stefanos
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
from itertools import chain

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui

from itaxotools.common.widgets import VLineSeparator, PushButton


class RecordFlag(Enum):
    Void = auto()
    Info = auto()
    Warn = auto()
    Fail = auto()


class RecordData:
    export_name: str = ''

    def __init__(self, data: object, formatter: Optional[str] = None):
        self.data = data
        if formatter:
            self.export_name = formatter.format(self.export_name)

    def export(self, path: Path) -> None:
        raise NotImplementedError()

    def view(self) -> QtWidgets.QWidget:
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

    def sorted(self):
        return list(chain(
            iter(record for record in self if record.type == RecordFlag.Fail),
            iter(record for record in self if record.type == RecordFlag.Warn),
            iter(record for record in self if record.type == RecordFlag.Info),
        ))

class RecordLogModel(QtCore.QAbstractListModel):

    def __init__(self, log):
        super().__init__()
        self.log = log

    def rowCount(self, parent=None):
        return len(self.log.records)

    def data(self, index, role):
        if (
            index.row() < 0 or
            index.row() >= len(self.log.records) or
            index.column() != 0
        ):
            return None
        record = self.log.sorted()[index.row()]
        if role == QtCore.Qt.DisplayRole:
            return str(record)
        if role == QtCore.Qt.UserRole:
            return record
        return None


class RecordLogDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._hovering = False

    def sizeHint(self, option, index):
        record = index.data(QtCore.Qt.UserRole)
        font_metrics = QtGui.QFontMetrics(option.font)
        width = font_metrics.horizontalAdvance(record.title)
        height = font_metrics.height()
        return QtCore.QSize(width + 28, height + 12)

    def paint(self, painter, option, index):
        self.initStyleOption(option, index)
        record = index.data(QtCore.Qt.UserRole)
        painter.save()
        deco = {
            RecordFlag.Void: '\u27A4',
            RecordFlag.Info: '\u2714',
            RecordFlag.Warn: '\u2718',
            RecordFlag.Fail: '\u2718',
        }[record.type]
        decoRect = QtCore.QRect(option.rect)
        decoRect.adjust(8, 0, 0, 0)
        painter.drawText(decoRect, QtCore.Qt.AlignVCenter, deco)
        if option.state & QtWidgets.QStyle.State_MouseOver and self._hovering:
            font = painter.font()
            font.setUnderline(True)
            painter.setFont(font)
        textRect = QtCore.QRect(option.rect)
        textRect.adjust(24, 0, 0, 0)
        painter.drawText(textRect, QtCore.Qt.AlignVCenter, record.title)
        painter.restore()

    def editorEvent(self, event, model, option, index):
        if event.type() == QtCore.QEvent.Type.MouseMove:
            hovering = (
                event.x() >= 24 and
                event.x() <= self.sizeHint(option, index).width() and
                event.y() <= self.parent().sizeHint().height())
            if hovering != self._hovering:
                if hovering:
                    self.parent().setCursor(QtCore.Qt.PointingHandCursor)
                else:
                    self.parent().unsetCursor()
                self._hovering = hovering
                # self.parent().update()
        elif (
            event.type() == QtCore.QEvent.MouseButtonRelease and
            event.button() == QtCore.Qt.LeftButton and
            self._hovering
        ):
            self.parent()._clicked(index)
        return super().event(event)


class RecordLogView(QtWidgets.QListView):
    link_clicked = QtCore.Signal(object)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMouseTracking(True)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Minimum,
            QtWidgets.QSizePolicy.Policy.Minimum)
        self.setStyleSheet("""
            RecordLogView {
                background: transparent;
                border: 0px solid transparent;
                outline: none;
                padding: 0px;
            }
        """)
        self.setSpacing(0)
        self.delegate = RecordLogDelegate(self)
        self.setItemDelegate(self.delegate)

    def sizeHint(self):
        width = self.sizeHintForColumn(0)
        height = self.sizeHintForRow(0)
        rows = 0
        if self.model():
            rows = self.model().rowCount()
        return QtCore.QSize(width, height * rows)

    def setLog(self, log):
        model = RecordLogModel(log)
        self.setModel(model)
        self.updateGeometry()

    def _clicked(self, index):
        record = index.data(QtCore.Qt.UserRole)
        self.link_clicked.emit(record)


class RecordDialog(QtWidgets.QDialog):
    def __init__(self, record: Record, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setWindowFlag(QtCore.Qt.WindowMaximizeButtonHint, True)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle(self.parent().title)
        self.record = record

        self.title = QtWidgets.QLabel(self.record.title)
        self.title.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        self.title.setStyleSheet("font-weight: bold;")

        self.description = QtWidgets.QLabel(self.record.description)
        self.description.setVisible(bool(record.description))
        self.description.setTextFormat(QtCore.Qt.RichText)
        self.description.setOpenExternalLinks(True)
        self.description.setTextInteractionFlags(
            QtCore.Qt.TextSelectableByMouse |
            QtCore.Qt.LinksAccessibleByMouse)
        self.description.setWordWrap(True)

        if self.record.data is not None:
            self.view = self.record.data.view()
        else:
            self.view = QtWidgets.QWidget()

        deco = {
            RecordFlag.Void: '\u27A4',
            RecordFlag.Info: '\u2714',
            RecordFlag.Warn: '\u2718',
            RecordFlag.Fail: '\u2718',
        }[self.record.type]
        self.deco = QtWidgets.QLabel(deco)

        self.ok = PushButton('OK')
        self.ok.clicked.connect(self.accept)
        self.ok.setDefault(True)
        self.save = PushButton('Save')
        self.save.clicked.connect(self.handleSave)
        self.save.setVisible(self.record.data is not None)

        buttons = QtWidgets.QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(self.save)
        buttons.addWidget(self.ok)
        buttons.setSpacing(8)
        buttons.setContentsMargins(0, 0, 0, 0)

        body = QtWidgets.QGridLayout()
        body.setRowMinimumHeight(0, 8)
        body.addWidget(self.deco, 1, 0, 1, 1)
        body.addWidget(self.title, 1, 1, 1, 2)
        body.setRowMinimumHeight(2, 8)
        body.addWidget(self.description, 3, 0, 1, 3)
        body.setRowMinimumHeight(4, 8)
        body.addWidget(self.view, 5, 0, 1, 3)
        body.setRowStretch(5, 1)
        body.setRowMinimumHeight(6, 24)
        body.setColumnStretch(2, 1)
        body.setHorizontalSpacing(6)
        body.setContentsMargins(0, 0, 0, 0)

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(body)
        layout.addLayout(buttons)
        layout.setContentsMargins(24, 16, 24, 16)
        self.setLayout(layout)

    def showEvent(self, event):
        super().showEvent(event)
        self.ok.setFocus()

    def handleSave(self):
        (path, _) = QtWidgets.QFileDialog.getSaveFileName(
            self, self.record.title + ' - Save',
            QtCore.QDir.currentPath() + '/' + self.record.data.export_name)
        if path:
            self.record.data.export(Path(path))
            self.accept()
