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

"""StepFilter"""

from PySide6 import QtCore
from PySide6 import QtWidgets

from lorem_text import lorem
from random import randint

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from .. import widgets
from .. import step_state_machine as ssm


class FilterItem(widgets.WidgetItem):
    map = {
        'name': 0,
        'filter': 1,
        'samples': 2,
        'nucleotides': 3,
        'uniform': 4,
        'missing': 5,
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.file = None
        alignment = QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter
        for col in range(1, 6):
            self.setTextAlignment(col, alignment)
        # font = self.font(0)
        # font.setBold(True)
        # self.setFont(0, font)
        self.name = lorem.words(randint(2, 6)).replace(' ', '_')
        self.filter = '-'
        self.samples = randint(6, 300)
        self.nucleotides = randint(200, 3000000)
        self.uniform = ['Yes', 'No'][randint(0, 1)]
        self.missing = randint(0, 9999) / 10000


class DataObject(object):
    pass


class StepFilter(ssm.StepState):

    title = 'Filter Character Sets'
    description = 'Rename and remove character sets'

    signalAdd = QtCore.Signal()
    signalDone = QtCore.Signal()
    signalRefresh = QtCore.Signal()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = DataObject()
        self.add_dummy_contents()

    def add_dummy_contents(self):
        count = randint(5, 50)
        for i in range(0, count):
            FilterItem(self.view)
        self.view.resizeColumnsToContents()
        self.sets.setValue(count)

    def work(self):
        time = randint(500, 2000)
        for i in range(0, int(time/10)):
            QtCore.QThread.msleep(10)
            self.worker.check()

    def draw(self):
        widget = QtWidgets.QWidget()

        text = ('Quisque tortor est, porttitor sed viverra ut, '
                'pharetra at nunc. Aenean vel congue dui. '
                'Vivamus auctor, quam se. \n'
                'Lorem ipsum dolor sit amet, consectetur adipiscing elit.'
                )
        label = QtWidgets.QLabel(text)

        frame = self.draw_frame()

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(label)
        layout.addWidget(frame, 1)
        layout.setSpacing(24)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)

        return widget

    def draw_summary(self):
        sets = widgets.InfoLabel('Total Sets')
        renamed = widgets.InfoLabel('Renamed')
        removed = widgets.InfoLabel('Removed')
        samples = widgets.InfoLabel('Unique Samples')

        for item in [sets, renamed, removed, samples]:
            item.setToolTip(lorem.words(randint(5, 15)))

        summary = QtWidgets.QHBoxLayout()
        summary.addWidget(sets)
        summary.addWidget(renamed)
        summary.addWidget(removed)
        summary.addWidget(samples)
        summary.addStretch(1)
        summary.setSpacing(24)
        summary.setContentsMargins(4, 0, 4, 0)

        self.sets = sets
        self.renamed = renamed
        self.removed = removed
        self.samples = samples

        return summary

    def draw_frame(self):
        view = widgets.TreeWidget()
        view.itemSelectionChanged.connect(self.handleItemSelectionChanged)
        view.setIndentation(0)
        view.setColumnCount(6)
        view.setHeaderLabels([
            ' Name', ' Filter', ' Samples',
            ' Nucleotides', ' Uniform', ' Missing'])

        headerItem = view.headerItem()
        headerItem.setToolTip(0, lorem.words(13))
        for col in range(1, 6):
            headerItem.setToolTip(col, lorem.words(randint(5, 15)))

        frame = common.widgets.Frame()

        rename = common.widgets.PushButton('Rename')
        rename.clicked.connect(self.handleRename)
        remove = common.widgets.PushButton('Remove')
        remove.clicked.connect(self.handleRemove)
        revert = common.widgets.PushButton('Revert')
        revert.clicked.connect(self.handleRevert)

        search = widgets.ViewSearchWidget(self, view)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(rename)
        controls.addWidget(remove)
        controls.addWidget(revert)
        controls.addStretch(1)
        controls.addWidget(search)
        controls.setSpacing(8)
        controls.setContentsMargins(0, 0, 0, 0)

        summary = self.draw_summary()

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(controls)
        layout.addWidget(view, 1)
        layout.addLayout(summary)
        layout.setSpacing(12)
        layout.setContentsMargins(16, 16, 16, 12)
        frame.setLayout(layout)

        self.view = view
        self.frame = frame
        self.remove = remove
        self.search = search

        return frame

    def refresh_contents(self):
        self.files.setValue(len(self.data.files))
        self.view.resizeColumnsToContents()

    def handleRename(self, checked=False):
        pass

    def handleRemove(self, checked=False):
        pass

    def handleRevert(self, checked=False):
        pass

    def handleItemSelectionChanged(self):
        pass
