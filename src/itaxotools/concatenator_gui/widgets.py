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

    def setValue(self, value='-'):
        self.value = value
        if isinstance(value, int):
            value = f'{value:,}'
        self.setText(f'{self.prefix}: {value}')


class HtmlLabel(QtWidgets.QLabel):
    def __init__(self, path):
        with open(path) as file:
            text = file.read()
        super().__init__(text)
        self.setTextFormat(QtCore.Qt.RichText)
        self.setOpenExternalLinks(True)
        self.setTextInteractionFlags(
            QtCore.Qt.TextSelectableByMouse |
            QtCore.Qt.LinksAccessibleByMouse)
        self.setWordWrap(True)


class _WidgetItem_meta(type(QtWidgets.QTreeWidgetItem)):
    def __new__(cls, name, bases, cdict):
        if 'fields' in cdict and isinstance(cdict['fields'], list):
            cdict['map'] = {v: i for i, v in enumerate(cdict['fields'])}
            cdict['unmap'] = {v: k for k, v in cdict['map'].items()}
        return super().__new__(cls, name, bases, cdict)


class WidgetItem(QtWidgets.QTreeWidgetItem, metaclass=_WidgetItem_meta):
    fields = []
    map = {}
    unmap = {}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.copyTextAlignment()

    def copyTextAlignment(self):
        if not self.treeWidget():
            return
        headerItem = self.treeWidget().headerItem()
        if not headerItem:
            return
        for column in range(0, headerItem.columnCount()):
            self.setTextAlignment(column, headerItem.textAlignment(column))

    def __setattr__(self, attr, value):
        super().__setattr__(attr, value)
        self.updateField(attr)

    def updateField(self, field):
        if field in self.fields:
            value = getattr(self, field)
            if isinstance(value, int):
                value = f'{value:,}'
            elif isinstance(value, float):
                value *= 100
                value = f'{value:.2f}%'
            elif isinstance(value, set):
                value = len(value)
            self.setText(self.map[field], str(value))
            if self.map[field] == 0:
                self.setToolTip(0, str(value))

    def __lt__(self, other):
        col = self.treeWidget().sortColumn()
        this = getattr(self, self.unmap[col])
        that = getattr(other, self.unmap[col])
        if isinstance(this, str) and isinstance(that, str):
            this = this.casefold()
            that = that.casefold()
        if isinstance(this, set) and isinstance(that, set):
            this = len(this)
            that = len(that)
        return this < that

    def setData(self, column, role, value):
        super().setData(column, role, value)
        if role == QtCore.Qt.EditRole and column in self.unmap:
            setattr(self, self.unmap[column], value)
            self.treeWidget().scrollToItem(self)


class ModelItem(WidgetItem):
    # Is it really needed to link model attibutes like this?
    def __init__(self, parent, model):
        super().__init__(parent)
        super().__setattr__('model', model)
        model_fields = [f for f in self.fields if hasattr(model, f)]
        for field in model_fields:
            self.updateField(field)

    def __setattr__(self, attr, value):
        super().__setattr__(attr, value)
        if hasattr(self.model, attr):
            setattr(self.model, attr, value)

    def __getattr__(self, attr):
        if hasattr(self.model, attr):
            return getattr(self.model, attr)
        return super().__getattribute__(attr)

    @property
    def samples_len(self):
        return len(self.model.samples)


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
            margin -= QtCore.QMargins(8, 0, 20, 0)

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

    def showEvent(self, event):
        item = self.topLevelItem(0)
        index = self.indexFromItem(item)
        self.setCurrentIndex(index)
        self.scrollToItem(item)

    def setColumnCount(self, columns, left_align=1):
        super().setColumnCount(columns)
        headerItem = self.headerItem()
        alignment = QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        for col in range(0, left_align):
            headerItem.setTextAlignment(col, alignment)
        alignment = QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter
        for col in range(left_align, columns):
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

    def iterate(self, start=None, end=None):
        assert self.topLevelItemCount() > 0
        if start is None:
            start = self.topLevelItem(0)
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
        self.setToolTip('Find next: F3')

        action = QtGui.QAction('Find next: F3', self)
        pixmap = common.widgets.VectorPixmap(
            common.resources.get('icons/svg/search.svg'),
            colormap=state.machine().parent().colormap_icon_light)
        action.setIcon(pixmap)
        action.setShortcut(QtGui.QKeySequence.FindNext)
        action.setStatusTip('Find next: F3')
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
