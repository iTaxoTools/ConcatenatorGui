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

"""StepDone"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui

from pathlib import Path

from itaxotools import common
import itaxotools.common.resources # noqa

from ..records import RecordLogView, RecordDialog
from ..diagnoser import SummaryReportView
from .. import step_state_machine as ssm
from .. import widgets

from . import wait


class StepDone(ssm.StepState):

    title = 'Results exported successfully'
    description = 'Concatenation complete'

    def cog(self):
        super().cog()
        machine = self.machine()
        transition = machine.navigateTransitionClear()
        transition.setTargetState(machine.states.input)
        self.addTransition(transition)
        self.transitions.new = transition

    def draw(self):
        widget = QtWidgets.QWidget()

        self.progress = wait.ProgressBar()
        self.progress.bar.setMaximum(1)
        self.progress.bar.setValue(1)

        self.confirm = self.progress.label
        self.confirm.setTextFormat(QtCore.Qt.RichText)
        self.confirm.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)

        self.diagnostics = QtWidgets.QLabel('Below you may find the summary reports for all performed diagnostics:')
        self.report_view = SummaryReportView()
        self.log_view = RecordLogView()

        self.report_view.clicked.connect(self.open_record)
        self.log_view.clicked.connect(self.open_record)

        path = common.resources.get(
            'itaxotools.concatenator_gui', 'docs/report.html')
        with open(path) as f:
            text = f.read()
        self.report = widgets.HtmlLabel(path)
        self.report.setOpenExternalLinks(False)
        self.report.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)
        self.report.linkActivated.connect(self.open_report_link)

        path = common.resources.get(
            'itaxotools.concatenator_gui', 'docs/done.html')
        self.label = widgets.HtmlLabel(path)
        self.label.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)
        self.label.setVisible(False)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.progress)
        layout.addWidget(self.diagnostics)
        layout.addSpacing(16)
        layout.addWidget(self.report_view)
        layout.addWidget(self.log_view)
        # layout.addWidget(self.report)
        # layout.addWidget(self.label)
        layout.addStretch(1)
        layout.setSpacing(16)
        layout.setContentsMargins(0, 0, 0, 32)
        widget.setLayout(layout)

        return widget

    def open_report_link(self, link):
        dir = self.machine().states.export.data.diagnoser.export_dir()
        file = Path(dir) / link
        url = QtCore.QUrl.fromLocalFile(str(file))
        QtGui.QDesktopServices.openUrl(url)

    def open_record(self, record):
        self.dialog = RecordDialog(record, self.machine().parent())
        self.dialog.setModal(True)
        self.dialog.show()

    def updateLabels(self):
        path = self.machine().states.export.data.target
        count_seqs = self.machine().states.export.data.seqs
        s = 's' if count_seqs > 1 else ''
        text = f'{count_seqs} sequence{s}'
        count_trees = self.machine().states.export.data.trees
        if count_trees:
            s = 's' if count_trees > 1 else ''
            text += f' and {count_trees} tree{s}'
        self.confirm.setText((
            f'<b>Successfully exported {text} to "{path.name}"</b>'))

    def warnDisjoint(self):
        if not self.machine().states.export.data.diagnoser:
            return
        group_count = self.machine().states.export.data.diagnoser.disjoint_groups
        if not group_count:
            return
        if group_count > 1:
            msgBox = QtWidgets.QMessageBox(self.machine().parent())
            msgBox.setWindowTitle(self.machine().parent().title)
            msgBox.setIcon(QtWidgets.QMessageBox.Warning)
            msgBox.setText(f'Disjoint sample groups detected!')
            msgBox.setInformativeText(
                f'We detected {group_count} distinct sample groups. '
                'This could be the result of slight differences in sample names between input files. '
                'Please open the disjoint group report to verify this is not a mistake.'
                )
            msgBox.setStandardButtons(
                QtWidgets.QMessageBox.Ignore | QtWidgets.QMessageBox.Open)
            msgBox.setDefaultButton(QtWidgets.QMessageBox.Ignore)
            button = self.machine().parent().msgShow(msgBox)
            if button == QtWidgets.QMessageBox.Open:
                self.open_report_link('disjoint_groups.txt')

    def onEntry(self, event):
        super().onEntry(event)
        self.updateLabels()
        # self.warnDisjoint()
        diagnoser = self.machine().states.export.data.diagnoser
        # for name, record in diagnoser.get_summary_report().records.items():
        #     print(str(record))
        #     print(record.data.data.to_string())
        #     print('')
        # print(diagnoser.get_record_log())
        self.report_view.setReport(diagnoser.get_summary_report())
        self.log_view.setLog(diagnoser.get_record_log())
