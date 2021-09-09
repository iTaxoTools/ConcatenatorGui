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

"""Custom widgets"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui

from time import time_ns

import re

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa


class SpinningCircle(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.handleTimer)
        self.timerStep = 10
        self.radius = 8
        self.period = 2
        self.span = 120
        self.width = 2

    def start(self):
        self.timer.start(self.timerStep)

    def stop(self):
        self.timer.stop()

    def handleTimer(self):
        self.repaint()

    def sizeHint(self):
        diameter = (self.radius + self.width) * 2
        return QtCore.QSize(diameter, diameter)

    def paintEvent(self, event):
        painter = QtGui.QPainter()
        painter.begin(self)

        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setBrush(QtCore.Qt.NoBrush)

        x = self.size().width()/2
        y = self.size().height()/2
        painter.translate(QtCore.QPoint(x, y))

        palette = QtGui.QGuiApplication.palette()
        weak = palette.color(QtGui.QPalette.Mid)
        bold = palette.color(QtGui.QPalette.Shadow)

        rad = self.radius
        rect = QtCore.QRect(-rad, -rad, 2 * rad, 2 * rad)

        painter.setPen(QtGui.QPen(weak, self.width, QtCore.Qt.SolidLine))
        painter.drawEllipse(rect)

        period_ns = int(self.period * 10**9)
        ns = time_ns() % period_ns
        degrees = - 360 * ns / period_ns
        painter.setPen(QtGui.QPen(bold, self.width, QtCore.Qt.SolidLine))
        painter.drawArc(rect, degrees * 16, self.span * 16)

        painter.end()


class InfoLabel(QtWidgets.QLabel):
    def __init__(self, text, value='-'):
        super().__init__()
        self.prefix = text
        self.setValue(value)

    def setValue(self, value):
        self.value = value
        if isinstance(value, int):
            value = f'{value:,}'
        self.setText(f'{self.prefix}: {value}')


class _WidgetItem_meta(type(QtWidgets.QTreeWidgetItem)):
    def __init__(cls, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if hasattr(cls, 'map') and isinstance(cls.map, dict):
            cls.unmap = {v: k for k, v in cls.map.items()}


class WidgetItem(QtWidgets.QTreeWidgetItem, metaclass=_WidgetItem_meta):
    map = {}

    def __setattr__(self, attr, value):
        super().__setattr__(attr, value)
        if attr in self.map:
            if isinstance(value, int):
                value = f'{value:,}'
            elif isinstance(value, float):
                value *= 100
                value = f'{value:.2f}%'
            self.setText(self.map[attr], str(value))
            self.setToolTip(self.map[attr], str(value))

    def __lt__(self, other):
        col = self.treeWidget().sortColumn()
        this = getattr(self, self.unmap[col])
        that = getattr(other, self.unmap[col])
        if isinstance(this, str) and isinstance(that, str):
            this = this.casefold()
            that = that.casefold()
        return this < that


class HeaderView(QtWidgets.QHeaderView):

    indicator_polygon = QtGui.QPolygon([
        QtCore.QPoint(-4, -1),
        QtCore.QPoint(0, 3),
        QtCore.QPoint(4, -1)])

    def paintSection(self, painter, rect, index):

        option = QtWidgets.QStyleOptionHeader()
        self.initStyleOptionForIndex(option, index)
        soh = QtWidgets.QStyleOptionHeader

        palette = QtGui.QGuiApplication.palette()
        dark = palette.color(QtGui.QPalette.Dark)
        window = palette.color(QtGui.QPalette.Window)
        light = palette.color(QtGui.QPalette.Light)
        midlight = palette.color(QtGui.QPalette.Midlight)
        mid = palette.color(QtGui.QPalette.Mid)
        black = palette.color(QtGui.QPalette.Text)

        top = rect.center() + QtCore.QPoint(0, - 2*rect.height())
        bot = rect.center() + QtCore.QPoint(0, rect.height())
        gradient = QtGui.QLinearGradient(top, bot)
        gradient.setColorAt(0, light)
        gradient.setColorAt(1, window)
        painter.fillRect(rect, gradient)

        painter.setPen(QtGui.QPen(midlight, 1, QtCore.Qt.SolidLine))
        painter.drawLine(rect.topRight(), rect.bottomRight())

        painter.setPen(QtGui.QPen(mid, 1, QtCore.Qt.SolidLine))
        painter.drawLine(rect.bottomLeft(), rect.bottomRight())

        margin = QtCore.QMargins(0, 0, 0, 0)
        if option.textAlignment & QtCore.Qt.AlignRight:
            margin -= QtCore.QMargins(20, 0, 8, 0)
        else:
            margin -= QtCore.QMargins(4, 0, 20, 0)
        # if option.position == soh.Beginning:
        # margin += QtCore.QMargins(4, 0, 0, 0)

        painter.setPen(QtGui.QPen(black, 1, QtCore.Qt.SolidLine))
        painter.drawText(rect + margin, option.textAlignment, option.text)

        if option.sortIndicator == soh.SortIndicator.None_:
            return

        painter.save()
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.translate(0, rect.center().y())
        if option.textAlignment & QtCore.Qt.AlignRight:
            painter.translate(rect.left() + 10, 0)
        else:
            painter.translate(rect.right() - 10, 0)
        if option.sortIndicator == soh.SortIndicator.SortDown:
            painter.scale(1, -1)
        painter.setBrush(QtCore.Qt.NoBrush)
        painter.setPen(QtGui.QPen(dark, 1.4, QtCore.Qt.SolidLine))
        painter.drawPolyline(self.indicator_polygon)
        painter.restore()

    def sectionSizeHint(self, column):
        option = QtWidgets.QStyleOptionHeader()
        self.initStyleOptionForIndex(option, column)
        m = self.fontMetrics().horizontalAdvance(option.text)
        return m + 30


class TreeWidget(QtWidgets.QTreeWidget):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setHeader(HeaderView(QtCore.Qt.Horizontal, self))
        self.setUniformRowHeights(True)
        self.header().setSectionsMovable(True)
        self.header().setStretchLastSection(False)
        self.header().setCascadingSectionResizes(False)
        self.header().setMinimumSectionSize(20)
        self.header().setSectionResizeMode(QtWidgets.QHeaderView.Fixed)
        self.header().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.setSelectionMode(QtWidgets.QTreeWidget.ExtendedSelection)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.setAllColumnsShowFocus(True)
        self.setAlternatingRowColors(True)
        self.setSortingEnabled(True)
        self.setStyleSheet("""
            QTreeView { border: 1px solid Palette(Mid); }
            QHeaderView::section {
                padding: 2px;
                }
            QTreeView::item {
                border: 0px;
                padding: 2px 6px;
                }
            QTreeView::item:selected {
                background: Palette(Highlight);
                color: Palette(Light);
                }
            """)

    def setColumnCount(self, columns):
        super().setColumnCount(columns)
        headerItem = self.headerItem()
        alignment = QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        headerItem.setTextAlignment(0, alignment)
        alignment = QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter
        for col in range(1, 6):
            headerItem.setTextAlignment(col, alignment)
        self.sortByColumn(0, QtCore.Qt.AscendingOrder)

    def check_item(self, item):
        if not isinstance(item, QtWidgets.QTreeWidgetItem):
            raise ValueError('Invalid item type.')
        if item.treeWidget() is not self:
            raise ValueError('Item does not belong to tree.')

    def get_next_item(self, item):
        self.check_item(item)
        if item.childCount() > 0:
            return item.child(0)
        parent = item.parent()
        while parent:
            index = parent.indexOfChild(item) + 1
            count = parent.childCount()
            if index < count:
                return parent.child(index)
            item = parent
            parent = item.parent()
        index = self.indexOfTopLevelItem(item) + 1
        count = self.topLevelItemCount()
        if index >= count:
            index = 0
        return self.topLevelItem(index)

    def iterate(self, start, end=None):
        if end is None:
            end = start
        self.check_item(start)
        self.check_item(end)
        yield start
        next_item = self.get_next_item(start)
        while next_item is not end:
            yield next_item
            next_item = self.get_next_item(next_item)
        yield end

    def resizeColumnToContents(self, column):
        if column < 0 or column >= self.header().count():
            return
        contents = self.sizeHintForColumn(column)
        header = 0
        if not self.header().isHidden():
            header = self.header().sectionSizeHint(column)
        self.header().resizeSection(column, max(contents, header))

    def resizeColumnsToContents(self):
        for column in range(1, self.header().count()):
            self.resizeColumnToContents(column)


class ViewSearchWidget(common.widgets.SearchWidget):

    def __init__(self, state, view, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.view = view
        self.state = state

        action = QtGui.QAction('Search', self)
        pixmap = common.widgets.VectorPixmap(
            common.resources.get_icon('search.svg'),
            colormap=state.machine().parent().colormap_icon_light)
        action.setIcon(pixmap)
        action.setShortcut(QtGui.QKeySequence.FindNext)
        action.setStatusTip('Search')
        action.triggered.connect(self.handleSearch)
        self.setSearchAction(action)

    def handleSearch(self, checked=False):
        found = None
        what = self.text()
        if not what:
            return

        current = self.view.currentItem()
        if current:
            next_item = self.view.get_next_item(current)
        else:
            current = self.view.topLevelItem(0)
            if not current:
                return
            next_item = current

        for item in self.view.iterate(next_item, current):
            text = item.text(0)
            if re.search(what, text, re.IGNORECASE):
                found = item
                break

        if found:
            self.view.setCurrentItem(found)
            self.view.scrollToItem(found)
