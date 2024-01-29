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

"""StepTrim"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui

from tempfile import TemporaryDirectory
from itertools import chain
from typing import Optional
from pathlib import Path
from time import sleep

from itaxotools import common
import itaxotools.common.widgets
import itaxotools.common.resources # noqa

from itaxotools.common.utility import AttrDict
from itaxotools.common.widgets import PushButton
from itaxotools.concatenator import (
    FileType, FileFormat, GeneStream, FileWriter,
    get_writer, get_extension, read_from_path)
from itaxotools.concatenator.library.model import GeneSeries, Operator
from itaxotools.concatenator.library.operators import OpUpdateMetadata
from itaxotools.concatenator.library.codons import (
    GeneticCode, ReadingFrame, _GC_DESCRIPTIONS)
from itaxotools.concatenator.library.operators import (
    OpChainGenes, OpTranslateGenes, OpApplyToGene, OpTagSet,
    OpUpdateMetadata, OpFilterGenes, OpGeneralInfo,
    OpGeneralInfoPerFile, OpGeneralInfoPerGene,
    OpGeneralInfoTagMafftRealigned,
    OpGeneralInfoTagPaddedLength,
    OpGeneralInfoTagPaddedCodonPosition)

from itaxotools.pygblocks import Options

from .. import model
from .. import widgets
from .. import step_state_machine as ssm
from ..trim import OpTrimGblocks, OpTrimClipKit, ClipKitTrimmingMode

from .wait import StepWaitBar
from .align import RichRadioButton


def expand_gap_characters(gc):
    return ''.join(c.upper() + c.lower() if c.isalpha() else c for c in gc)


class TrimItem(widgets.ModelItem):
    fields = [
        'display_name',
        'action',
        'samples_len',
        'nucleotides',
        'missing',
        'uniform',
        'aligned_string',
        ]

    def __init__(self, parent, charset: model.Charset):
        super().__init__(parent, charset)
        self.setFlags(QtCore.Qt.ItemIsSelectable |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemNeverHasChildren)
        self.samples_len = len(self.model.samples)
        self.tag_new('trimmed')
        self.refresh()

    def trim(self):
        self.trimmed = True
        self.refresh()

    def clear(self):
        self.trimmed = False
        self.refresh()

    def toggle(self):
        self.trimmed = not self.trimmed
        self.refresh()

    def refresh(self):
        self.updateField('action')
        self.updateField('aligned_string')
        self.tag_set('trimmed', self.trimmed)
        self.treeWidget().signalTagUpdate.emit()
        self.setBold(self.trimmed)

    def setBold(self, value):
        font = self.font(0)
        font.setBold(value)
        self.setFont(0, font)

    @property
    def action(self):
        return 'Trim' if self.trimmed else '-'

    @property
    def aligned_string(self):
        return 'Yes' if self.aligned else 'No'


class GblocksOptions(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._option_count = -1

        options = QtWidgets.QGridLayout()
        options.setHorizontalSpacing(16)
        options.setVerticalSpacing(12)
        options.setColumnStretch(99, 1)

        self.fields = AttrDict()

        self.fields.GT = self.getPercentageField()
        self.addOption(options, self.fields.GT, 'Gaps allowed', 'Maximum percentage of gaps allowed for any position.')

        self.fields.IS = self.getPercentageField()
        self.addOption(options, self.fields.IS, 'Conservation threshold (IS)', 'Minimum percentage of identical sequences for a conserved position.')

        self.fields.FS = self.getPercentageField()
        self.addOption(options, self.fields.FS, 'Flank position threshold (FS)', 'Minimum percentage of identical sequences for a flank position.')

        self.fields.GC = self.getTextField()
        self.fields.GC.setPlaceholderText('-?*')
        self.addOption(options, self.fields.GC, 'Gap definition', 'Definition of gap characters, such as -?*XxNn.')

        self.fields.CP = self.getIntegerField()
        self.fields.CP.setPlaceholderText('8')
        self.addOption(options, self.fields.CP, 'Nonconserved (CP)', 'Maximum number of contiguous nonconserved positions.')

        self.fields.BL = self.getIntegerField()
        self.fields.BL.setPlaceholderText('15')
        self.addOption(options, self.fields.BL, 'Block length (BL)', 'Minimum length of a block, for both 1st and 2nd iteration.')

        config_restore = PushButton('Restore defaults', onclick=self.restoreDefaults)
        config_restore.setFixedWidth(200)

        config_restore_layout = QtWidgets.QHBoxLayout()
        config_restore_layout.addWidget(config_restore, 1)
        config_restore_layout.addStretch(1)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 0, 0, 0)
        layout.setSpacing(28)
        layout.addLayout(options)
        layout.addLayout(config_restore_layout)

        self.restoreDefaults()

    def addOption(self, layout: QtWidgets.QGridLayout, field: QtWidgets.QWidget, name: str, tooltip: str):
        label = QtWidgets.QLabel(name + ':')
        label.setToolTip(tooltip)
        field.setToolTip(tooltip)

        self._option_count += 1

        col, row = divmod(self._option_count, 3)
        layout.addWidget(label, row, col * 3)
        layout.addWidget(field, row, col * 3 + 1)

        if row == 0:
            layout.setColumnMinimumWidth(col * 3 + 2, 48)

    def getPercentageField(self):
        spinbox = QtWidgets.QDoubleSpinBox()
        spinbox.setSuffix(' %')
        spinbox.setMinimum(0)
        spinbox.setMaximum(100)
        spinbox.setDecimals(2)
        spinbox.setSingleStep(1.0)
        spinbox.setFixedWidth(100)
        return spinbox

    def getIntegerField(self):
        edit = QtWidgets.QLineEdit()
        validator = QtGui.QIntValidator()
        edit.setValidator(validator)
        edit.setFixedWidth(100)
        return edit

    def getTextField(self):
        edit = QtWidgets.QLineEdit()
        edit.setFixedWidth(100)
        return edit

    def restoreDefaults(self):
        self.fields.IS.setValue(50.00)
        self.fields.FS.setValue(85.00)
        self.fields.GT.setValue(0.00)
        self.fields.CP.setText('8')
        self.fields.BL.setText('15')
        self.fields.GC.setText('-?*')

    def getTextOrPlaceholder(self, field):
        if field.text():
            return field.text()
        return field.placeholderText()

    def toOptions(self) -> Options:
        IS = self.fields.IS.value() / 100
        FS = self.fields.FS.value() / 100
        GT = self.fields.GT.value() / 100

        CP = int(self.getTextOrPlaceholder(self.fields.CP))
        BL = int(self.getTextOrPlaceholder(self.fields.BL))
        GC = self.getTextOrPlaceholder(self.fields.GC)

        return Options(
            IS=None,
            FS=None,
            GT=None,

            CP=CP,
            BL1=BL,
            BL2=BL,
            GC=expand_gap_characters(GC),

            IS_percent=IS,
            FS_percent=FS,
            GT_percent=GT,
        )


class ClipkitOptions(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._gappy_default = 90.00
        self._option_count = -1

        options = QtWidgets.QGridLayout()
        options.setHorizontalSpacing(16)
        options.setVerticalSpacing(12)
        options.setColumnStretch(99, 1)

        self.fields = AttrDict()

        self.fields.mode = self.getModeField()
        self.addOption(options, self.fields.mode, 'ClipKIT mode', 'One of the various trimming modes implemented in ClipKIT.', large=True)

        self.fields.gaps = self.getGappyField()
        self.addOption(options, self.fields.gaps, 'Gaps allowed', 'Maximum percentage of gaps allowed for any position.')

        self.fields.gap_characters = self.getTextField()
        self.fields.gap_characters.setPlaceholderText('-?*')
        self.addOption(options, self.fields.gap_characters, 'Gap definition', 'Definition of gap characters, such as -?*XxNn.')

        self.fields.mode.currentIndexChanged.connect(self.update_gappy_field)

        config_restore = PushButton('Restore defaults', onclick=self.restoreDefaults)
        config_restore.setFixedWidth(200)

        config_restore_layout = QtWidgets.QHBoxLayout()
        config_restore_layout.addWidget(config_restore, 1)
        config_restore_layout.addStretch(1)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 0, 0, 0)
        layout.setSpacing(28)
        layout.addLayout(options)
        layout.addLayout(config_restore_layout)

        self.restoreDefaults()

    def addOption(self, layout: QtWidgets.QGridLayout, field: QtWidgets.QWidget, name: str, tooltip: str, large=False):
        label = QtWidgets.QLabel(name + ':')
        label.setToolTip(tooltip)
        field.setToolTip(tooltip)

        self._option_count += 1
        row, col = divmod(self._option_count, 2)

        span = 1
        if large:
            self._option_count += 1
            span = 4

        layout.addWidget(label, row, col * 3, 1, 1)
        layout.addWidget(field, row, col * 3 + 1, 1, span)

        if row == 0:
            layout.setColumnMinimumWidth(col * 3 + 2, 48)

    def getModeField(self):
        combobox = QtWidgets.QComboBox()
        for mode in ClipKitTrimmingMode:
            combobox.addItem(f'{mode.mode}: {mode.description}', mode)
        return combobox

    def getGappyField(self):
        spinbox = QtWidgets.QDoubleSpinBox()
        spinbox.setSuffix(' %')
        spinbox.setMinimum(0)
        spinbox.setMaximum(100)
        spinbox.setDecimals(2)
        spinbox.setSingleStep(1.0)
        spinbox.setFixedWidth(100)
        return spinbox

    def getTextField(self):
        edit = QtWidgets.QLineEdit()
        edit.setFixedWidth(100)
        return edit

    def restoreDefaults(self):
        self.fields.mode.setCurrentIndex(0)
        self.fields.gaps.setValue(self._gappy_default)
        self.fields.gap_characters.setText('-?*')
        self.update_gappy_field(0)

    def update_gappy_field(self, index: int):
        mode = self.fields.mode.itemData(index)
        if not 'gappy' in mode.mode:
            self.fields.gaps.setMinimum(0)
            self.fields.gaps.setMaximum(0)
            self.fields.gaps.setValue(0)
            self.fields.gaps.setSpecialValueText('N/A')
            self.fields.gaps.setEnabled(False)
        else:
            self.fields.gaps.setMinimum(0)
            self.fields.gaps.setMaximum(100)
            self.fields.gaps.setValue(self._gappy_default)
            self.fields.gaps.setSpecialValueText('')
            self.fields.gaps.setEnabled(True)

    def getTextOrPlaceholder(self, field):
        if field.text():
            return field.text()
        return field.placeholderText()

    def toOptions(self) -> Options:
        mode = self.fields.mode.currentData()
        gaps = self.fields.gaps.value() / 100
        gap_characters = self.getTextOrPlaceholder(self.fields.gap_characters)

        options = dict(
            mode=mode.name,
            gap_characters=expand_gap_characters(gap_characters),
        )

        if 'gappy' in mode.name:
            options['gaps'] = gaps

        return options


class StepTrimOptions(ssm.StepState):

    title = 'Trim Sequence Sites'
    description = 'Select trimming toolkit'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data = AttrDict()
        self.data.skip = False
        self.data.last_method = None
        self.data.last_options = None

    def draw(self):
        widget = QtWidgets.QWidget()

        path = common.resources.get(
            'itaxotools.concatenator_gui', 'docs/trim.html')
        head_label = widgets.HtmlLabel(path)

        self.radios = AttrDict()
        self.radios.gblocks = RichRadioButton('pyGblocks:', 'eliminate poorly aligned positions and divergent regions.', widget)
        self.radios.clipkit = RichRadioButton('ClipKIT:', 'only keep phylogenetically informative sites.', widget)
        self.radios.skip = QtWidgets.QRadioButton('Skip trimming', widget)
        self.radios.skip.setStyleSheet("QRadioButton { letter-spacing: 1px; }")

        radios = QtWidgets.QVBoxLayout()
        radios.addWidget(self.radios.gblocks)
        radios.addWidget(self.radios.clipkit)
        radios.addSpacing(16)
        radios.addWidget(self.radios.skip)
        radios.setSpacing(16)
        radios.setContentsMargins(16, 0, 0, 0)

        self.options = AttrDict()
        self.options.gblocks = GblocksOptions()
        self.options.clipkit = ClipkitOptions()
        self.options.skip = QtWidgets.QWidget()

        self.options_label = QtWidgets.QLabel('You may configure the trimming parameters below. Hover options for more information.')

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(head_label)
        layout.addLayout(radios)
        layout.addSpacing(12)
        layout.addWidget(self.options_label)
        layout.addWidget(self.options.gblocks)
        layout.addWidget(self.options.clipkit)
        layout.addWidget(self.options.skip)
        layout.addStretch(1)
        layout.setSpacing(24)
        layout.setContentsMargins(0, 0, 0, 12)
        widget.setLayout(layout)

        for radio in self.radios:
            radio.toggled.connect(self.update_options_shown)

        self.radios.skip.setChecked(True)

        return widget

    def get_method(self):
        for method in self.radios.keys():
            if self.radios[method].isChecked():
                return method

    def get_options_dict(self):
        method = self.get_method()
        if method == 'gblocks':
            return self.options.gblocks.toOptions().as_dict()
        elif method == 'clipkit':
            return self.options.clipkit.toOptions()
        return {}

    def update_options_shown(self, status):
        if not status:
            return
        method = self.get_method()
        for option in self.options:
            option.setVisible(False)
        self.options[method].setVisible(True)
        self.options_label.setVisible(method != "skip")

    def onExit(self, event):
        super().onExit(event)
        self.data.skip = self.radios.skip.isChecked()
        method = self.get_method()
        options = self.get_options_dict()

        if self.data.last_method != method or self.data.last_options != options:
            self.data.last_method = method
            self.data.last_options = options
            self.timestamp_set()


class StepTrimSetsEdit(ssm.StepTriStateEdit):

    description = 'Select which markers to trim'

    def onEntry(self, event):
        super().onEntry(event)
        last_filter_update = self.machine().states.filter.timestamp_get()
        if last_filter_update > self.timestamp_get():
            self.populate_view()
            self.timestamp_set()
        for item in self.view.iterate():
            item.refresh()

    def populate_view(self):
        self.view.clear()
        charsets = [
            cs for cs in self.machine().states.input.data.charsets.values()
            if cs.translation is not None]
        for charset in charsets:
            TrimItem(self.view, charset)
        self.view.resizeColumnsToContents()
        self.sets.setValue(len(charsets))
        self.updateSummary()

    def draw(self):
        widget = QtWidgets.QWidget()

        text = (
            'Select markers for trimming by double-clicking, or by '
            'highlighting them and then clicking "Trim".')
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
        sets = widgets.InfoLabel('Markers')
        marked = widgets.InfoLabel('Selected', 0)

        sets.setToolTip('Total number of markers.')
        marked.setToolTip('Number of markers pending trimming.')

        summary = QtWidgets.QHBoxLayout()
        summary.addWidget(sets)
        summary.addWidget(marked)
        summary.addStretch(1)
        summary.setSpacing(24)
        summary.setContentsMargins(4, 0, 4, 0)

        self.sets = sets
        self.marked = marked

        return summary

    def draw_frame(self):
        frame = common.widgets.Frame()

        view = widgets.TreeWidget()
        view.signalTagUpdate.connect(self.updateSummary)
        view.itemActivated.connect(self.handleActivated)
        view.setIndentation(0)
        view.setColumnCount(7, 2)
        view.setHeaderLabels([
            'Name', 'Action', 'Samples',
            'Nucleotides', 'Missing', 'Uniform', 'Aligned'])

        headerItem = view.headerItem()
        headerItem.setToolTip(0, 'Marker name')
        headerItem.setToolTip(1, 'Pending action')
        headerItem.setToolTip(2, 'Total number of sequences')
        headerItem.setToolTip(3, 'Total number of nucleotide characters')
        headerItem.setToolTip(4, 'Proportion of missing data')
        headerItem.setToolTip(5, 'Are all sequences of the same length?')
        headerItem.setToolTip(6, 'Were the sequences aligned with MAFFT?')

        all = common.widgets.PushButton('Trim All', onclick=self.handleAll)
        trim = common.widgets.PushButton('Trim', onclick=self.handleTrim)
        clear = common.widgets.PushButton('Clear', onclick=self.handleClear)

        search = widgets.ViewSearchWidget(self, view)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(all)
        controls.addWidget(trim)
        controls.addWidget(clear)
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
        self.search = search

        return frame

    def updateSummary(self):
        self.marked.setValue(self.view.tag_get('trimmed'))

    def handleTrim(self, checked=False):
        item = None
        for item in self.view.selectedItems():
            item.trim()
        self.view.scrollToItem(item)

    def handleClear(self, checked=False):
        item = None
        for item in self.view.selectedItems():
            item.clear()
        self.view.scrollToItem(item)

    def handleAll(self, checked=False):
        self.view.selectAll()
        self.handleTrim()

    def handleActivated(self, item, column):
        item.toggle()
        self.view.scrollToItem(item)


class StepTrimSetsWait(StepWaitBar):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger.setWordWrapMode(QtGui.QTextOption.WrapAnywhere)


class StepTrimSetsDone(ssm.StepTriStateDone):
    description = 'Trimming complete'

    def onEntry(self, event):
        super().onEntry(event)
        if self.result:
            s = 's' if self.result > 1 else ''
            self.parent().update(
                text=f'Successfully trimmed {str(self.result)} marker{s}.')
        else:
            self.parent().update(
                text='Successfully trimmed sequences.')


class StepTrimSetsFail(ssm.StepTriStateFail):
    description = 'Trimming failed'

    def onEntry(self, event):
        super().onEntry(event)
        message = f'{type(self.exception).__name__}'
        self.parent().update(
            text=f'Position trimming failed: {message}')


class StepTrimSets(ssm.StepTriState):
    title = 'Trim Sequence Sites'

    StepEdit = StepTrimSetsEdit
    StepWait = StepTrimSetsWait
    StepDone = StepTrimSetsDone
    StepFail = StepTrimSetsFail

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.temp_cache = None
        self.charsets_cached = set()
        self.charsets_last = set()

    def onEntry(self, event):
        super().onEntry(event)
        input_update = self.machine().states.input.timestamp_get()
        last_filter_update = self.machine().states.filter.timestamp_get()
        align_sets_update = self.machine().states.align_sets.timestamp_get()
        trim_options_update = self.machine().states.trim_options.timestamp_get()
        last_update = max(input_update, last_filter_update, align_sets_update, trim_options_update)
        if last_update > self.timestamp_get():
            self.temp_cache = TemporaryDirectory(prefix='concat_trim_cache_')
            self.charsets_cached = set()

    def onExit(self, event):
        super().onExit(event)
        charsets = {
            k for k, v in self.machine().states.input.data.charsets.items()
            if v.trimmed and v.translation is not None}
        if self.charsets_last != charsets:
            self.charsets_last = charsets
            self.timestamp_set()

    def work(self):
        charsets = {
            k for k, v in self.machine().states.input.data.charsets.items()
            if v.trimmed and v.translation is not None and v.samples
            and v.name not in self.charsets_cached}
        method = self.machine().states.trim_options.get_method()
        with self.states.wait.redirect():
            stream = self.work_get_stream(charsets)
            self.work_trim(stream, method, len(charsets))
        return len(charsets)

    def work_nothing(self, method: str, total:int):
        if method == 'gblocks':
            title = 'pyGblocks'
        elif method == 'clipkit':
            title = 'ClipKIT'
        else:
            raise Exception('Unexpected trimming toolkit')

        print(f'Trimming sequences using {title}...')
        print()
        print(f'Found {total} markers already cached, skipping!')

    def work_trim(self, stream: GeneStream, method: str, total: int):
        if method == 'gblocks':
            title = 'pyGblocks'
            options = self.machine().states.trim_options.options.gblocks.toOptions()
            options_dict = options.as_dict()
            operator = OpTrimGblocks(options)
        elif method == 'clipkit':
            title = 'ClipKIT'
            options_dict = self.machine().states.trim_options.options.clipkit.toOptions()
            operator = OpTrimClipKit(options_dict)
        else:
            raise Exception('Unexpected trimming toolkit')

        self.update(0, 0, 'Getting ready...')
        print(f'Starting trimming for {total} markers '
              f'({len(self.charsets_cached)} already cached)...')
        print()
        print(f'Selected toolkit: {title}')
        print()
        print('Options:\n')
        for key, value in options_dict.items():
            print(f'- {key}: {value}')
        print(f'\n{"-"*20}\n')

        def cache_func(series):
            if series:
                self.charsets_cached.add(series.name)
            return series

        cache_operator = OpApplyToGene(cache_func)

        stream = stream.pipe(operator).pipe(cache_operator)

        path = Path(self.temp_cache.name)
        writer = get_writer(FileType.Directory, FileFormat.Fasta)
        writer.params.translate_missing.value = ''
        writer.params.translate_gap.value = ''
        writer.params.padding.value = ''
        writer.params.sanitize.value = False
        writer.params.drop_empty.value = False
        writer(stream, path,)

        print('Completed site trimming!')
        print()
        self.update(1, 1, 'text')

        # if total > 0:
        #     self.timestamp_set()

    def work_get_stream(self, charsets):
        input_streams = [
            read_from_path(file.path)
            for file in self.machine().states.input.data.files.values()]
        aligned_cache = Path(self.machine().states.align_sets.temp_cache.name)
        aligned_charsets = {
            k for k, v in self.machine().states.input.data.charsets.items()
            if v.aligned and v.translation is not None}
        aligned_stream = (
            read_from_path(aligned_cache)
            .pipe(OpFilterGenes(aligned_charsets))
            )
        all_streams = [aligned_stream] + input_streams
        stream = (
            GeneStream(chain(*all_streams))
            .pipe(OpChainGenes())
            .pipe(OpFilterGenes(charsets))
            )
        return stream

    def skipAll(self):
        skip = self.machine().states['trim_options'].data.skip
        return bool(skip)

    def skipWait(self):
        skip = self.states['edit'].marked.value == 0
        return bool(skip)

    def onCancel(self, exception):
        # self.process.quit()
        pass

    def onFail(self, exception, trace):
        # raise exception
        self.states.wait.logio.writeline('')
        self.states.wait.logio.writeline(trace)
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Critical)
        msgBox.setText(type(exception).__name__)
        msgBox.setInformativeText(str(exception))
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.machine().parent().msgShow(msgBox)

    def filterCancel(self, event):
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Question)
        msgBox.setText('Cancel trimming?')
        msgBox.setStandardButtons(
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.No)
        res = self.machine().parent().msgShow(msgBox)
        return res == QtWidgets.QMessageBox.Yes

    def filterNext(self, event):
        bad_charsets = {
            k for k, v in self.machine().states.input.data.charsets.items()
            if v.trimmed and v.uniform == 'No' and not v.aligned}
        if not bad_charsets:
            return True
        msgBox = QtWidgets.QMessageBox(self.machine().parent())
        msgBox.setWindowTitle(self.machine().parent().title)
        msgBox.setIcon(QtWidgets.QMessageBox.Critical)
        msgBox.setText('Cannot trim markers that are not \nof uniform length or aligned.')
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.machine().parent().msgShow(msgBox)
        return False
