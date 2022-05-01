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

from ..records import RecordLogView, RecordDialog, RecordFlag
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

        self.diagnostics = QtWidgets.QLabel('Below you may find the results of data validation:')
        self.report_view = SummaryReportView()
        self.log_view = RecordLogView()

        self.report_view.clicked.connect(self.open_record)
        self.log_view.clicked.connect(self.open_record)

        path = common.resources.get(
            'itaxotools.concatenator_gui', 'docs/done.html')
        self.label = widgets.HtmlLabel(path)
        self.label.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)

        views = QtWidgets.QVBoxLayout()
        views.addWidget(self.report_view)
        views.addWidget(self.log_view)
        views.setSpacing(16)
        views.setContentsMargins(16, 4, 16, 4)

        self.validation = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.diagnostics)
        layout.addSpacing(8)
        layout.addLayout(views)
        layout.setContentsMargins(0, 0, 0, 0)
        self.validation.setLayout(layout)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.progress)
        layout.addWidget(self.validation)
        layout.addWidget(self.label)
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
        text = f'{count_seqs} marker{s}'
        count_trees = self.machine().states.export.data.trees
        if count_trees:
            s = 's' if count_trees > 1 else ''
            text += f' and {count_trees} tree{s}'
        self.confirm.setText((
            f'<b>Successfully exported {text} to "{path.name}"</b>'))

    def notifyWarnings(self, log):
        warnings = (record for record in log if record.type != RecordFlag.Info)
        if not any(warnings):
            return

        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Warning)
        msgBox.setText('We detected some possible problems with the data set.')
        msgBox.setInformativeText('Please click on the issued warnings for details.')
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        button = self.machine().parent().msgShow(msgBox)

    def onEntry(self, event):
        super().onEntry(event)
        self.updateLabels()
        diagnoser = self.machine().states.export.data.diagnoser
        report = diagnoser.get_summary_report()
        record_log = diagnoser.get_record_log()
        self.report_view.setReport(report)
        self.log_view.setLog(record_log)
        self.notifyWarnings(record_log)

        showValidation = bool(any(record_log) or report is not None)
        self.validation.setVisible(showValidation)
