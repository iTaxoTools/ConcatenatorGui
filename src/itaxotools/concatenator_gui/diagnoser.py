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

"""Data diagnoser"""


from typing import Dict, List, Optional, Tuple, Union
from tempfile import TemporaryDirectory
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import json

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui

from itaxotools.common import resources
from itaxotools.common.utility import AttrDict
from itaxotools.concatenator.library.operators import (
    OpChainGenes, OpTranslateGenes, OpApplyToGene, OpTagSet,
    OpUpdateMetadata, OpFilterGenes, OpGeneralInfo,
    OpGeneralInfoPerFile, OpGeneralInfoPerGene,
    OpGeneralInfoTagMafftRealigned,
    OpGeneralInfoTagPaddedLength,
    OpGeneralInfoTagPaddedCodonPosition)

from .bouncer import OpSequenceBouncer
from .records import Record, RecordFlag, RecordLog, RecordData


@dataclass
class DiagnoserParams:
    report: bool = True
    disjoint: bool = True
    foreign: bool = True
    outliers: bool = False
    iqr: float = 20.0


class SeriesModel(QtCore.QAbstractTableModel):

    def __init__(self, series: pd.Series, parent=None):
        super().__init__(parent)
        self._series = series

    def rowCount(self, parent=QtCore.QModelIndex()) -> int:
        if parent == QtCore.QModelIndex():
            return len(self._series)
        return 0

    def columnCount(self, parent=QtCore.QModelIndex()) -> int:
        if parent == QtCore.QModelIndex():
            return 1
        return 0

    def data(self, index: QtCore.QModelIndex, role=QtCore.Qt.ItemDataRole):
        if not index.isValid():
            return None

        if role == QtCore.Qt.DisplayRole:
            value = self._series.iloc[index.row()]
            if isinstance(value, float):
                return f'{value:.2f}'
            return str(value)

        return None

    def headerData(
        self, section: int, orientation: QtCore.Qt.Orientation, role: QtCore.Qt.ItemDataRole
    ):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                return str(self._series.name)

            if orientation == QtCore.Qt.Vertical:
                return str(self._series.index[section])

        return None


class DataFrameModel(QtCore.QAbstractTableModel):

    def __init__(self, dataframe: pd.DataFrame, parent=None):
        super().__init__(parent)
        self._dataframe = dataframe

    def rowCount(self, parent=QtCore.QModelIndex()) -> int:
        if parent == QtCore.QModelIndex():
            return len(self._dataframe)

        return 0

    def columnCount(self, parent=QtCore.QModelIndex()) -> int:
        if parent == QtCore.QModelIndex():
            return len(self._dataframe.columns)

        return 0

    def data(self, index: QtCore.QModelIndex, role=QtCore.Qt.ItemDataRole):
        if not index.isValid():
            return None

        if role == QtCore.Qt.DisplayRole:
            value = self._dataframe.iloc[index.row(), index.column()]
            if isinstance(value, float):
                return f'{value:.2f}'
            return str(value)

        return None

    def headerData(
        self, section: int, orientation: QtCore.Qt.Orientation, role: QtCore.Qt.ItemDataRole
    ):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                return str(self._dataframe.columns[section])

            if orientation == QtCore.Qt.Vertical:
                return str(self._dataframe.index[section])

        return None


class ForeignPairsModel(QtCore.QAbstractTableModel):
    def __init__(self, pairs: List[Tuple[str, str]], parent=None):
        super().__init__(parent)
        self._pairs = pairs

    def rowCount(self, parent=QtCore.QModelIndex()) -> int:
        if parent == QtCore.QModelIndex():
            return len(self._pairs)

        return 0

    def columnCount(self, parent=QtCore.QModelIndex()) -> int:
        if parent == QtCore.QModelIndex():
            return 2

        return 0

    def data(self, index: QtCore.QModelIndex, role=QtCore.Qt.ItemDataRole):
        if not index.isValid():
            return None

        if role == QtCore.Qt.DisplayRole:
            pair = self._pairs[index.row()]
            return str(pair[index.column()])

        return None

    def headerData(
        self, section: int, orientation: QtCore.Qt.Orientation, role: QtCore.Qt.ItemDataRole
    ):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                return f'Sample {section + 1}'

            if orientation == QtCore.Qt.Vertical:
                return str(section + 1)

        return None


class DisjointGroupsModel(QtCore.QAbstractItemModel):
    def __init__(self, groups: List[List[str]], parent=None):
        super().__init__(parent)
        self._groups = groups

    def index(self, row: int, column: int, parent=QtCore.QModelIndex()) -> QtCore.QModelIndex:
        if not self.hasIndex(row, column, parent):
            return QtCore.QModelIndex()

        if column != 0:
            return QtCore.QModelIndex()

        if not parent.isValid():
            return self.createIndex(row, 0, self._groups)
        elif parent.internalPointer() is self._groups:
            return self.createIndex(row, 0, self._groups[parent.row()])
        return QtCore.QModelIndex()

    def parent(self, index=QtCore.QModelIndex()) -> QtCore.QModelIndex:
        if not index.isValid():
            return QtCore.QModelIndex()
        ptr = index.internalPointer()
        if ptr is self._groups:
            return QtCore.QModelIndex()
        try:
            pos = self._groups.index(ptr)
            return self.createIndex(pos, 0, self._groups)
        except ValueError:
            return QtCore.QModelIndex()
        return QtCore.QModelIndex()

    def rowCount(self, parent=QtCore.QModelIndex()) -> int:
        if not parent.isValid():
            return len(self._groups)
        elif parent.internalPointer() is self._groups:
            return len(self._groups[parent.row()])
        return 0

    def columnCount(self, parent=QtCore.QModelIndex()) -> int:
        return 1

    def data(self, index: QtCore.QModelIndex, role=QtCore.Qt.ItemDataRole):
        if not index.isValid():
            return None

        if role == QtCore.Qt.DisplayRole:
            ptr = index.internalPointer()
            if ptr is self._groups:
                return str(f'Group {index.row() + 1}')
            return str(ptr[index.row()])

        return None


class OutliersModel(QtCore.QAbstractItemModel):
    def __init__(self, outliers: Dict[str, Union[str, List[str]]], parent=None):
        super().__init__(parent)
        outliers = {k: self._normalize(v) for k, v in outliers.items() if v}
        self._genes = list(outliers.keys())
        self._outliers = list(outliers.values())

    def index(self, row: int, column: int, parent=QtCore.QModelIndex()) -> QtCore.QModelIndex:
        if not self.hasIndex(row, column, parent):
            return QtCore.QModelIndex()

        if column != 0:
            return QtCore.QModelIndex()

        if not parent.isValid():
            return self.createIndex(row, 0, self._outliers)
        elif parent.internalPointer() is self._outliers:
            return self.createIndex(row, 0, self._outliers[parent.row()])
        return QtCore.QModelIndex()

    def parent(self, index=QtCore.QModelIndex()) -> QtCore.QModelIndex:
        if not index.isValid():
            return QtCore.QModelIndex()
        ptr = index.internalPointer()
        if ptr is self._outliers:
            return QtCore.QModelIndex()
        try:
            pos = self._outliers.index(ptr)
            return self.createIndex(pos, 0, self._outliers)
        except ValueError:
            return QtCore.QModelIndex()
        return QtCore.QModelIndex()

    def rowCount(self, parent=QtCore.QModelIndex()) -> int:
        if not parent.isValid():
            return len(self._outliers)
        elif parent.internalPointer() is self._outliers:
            return len(self._outliers[parent.row()])
        return 0

    def columnCount(self, parent=QtCore.QModelIndex()) -> int:
        return 1

    def data(self, index: QtCore.QModelIndex, role=QtCore.Qt.ItemDataRole):
        if not index.isValid():
            return None

        if role == QtCore.Qt.DisplayRole:
            ptr = index.internalPointer()
            if ptr is self._outliers:
                return str(self._genes[index.row()])
            return str(ptr[index.row()])

        return None

    def _normalize(self, data: Union[str, List[str]]) -> List[str]:
        if isinstance(data, str):
            return [data]
        return data


class TableView(QtWidgets.QTableView):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.setStyleSheet("QHeaderView::section {padding: 2px 5px 2px 5px;}")

    def setModel(self, model):
        super().setModel(model)
        self.resizeColumnsToContents()
        self.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Fixed)
        for column in range(0, model.columnCount()):
            self.setColumnWidth(column, self.columnWidth(column) + 10)


class TreeView(QtWidgets.QTreeView):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)

    def setModel(self, model):
        super().setModel(model)
        self.header().hide()
        self.expandAll()


class RecordTable(RecordData):
    def _export_table(
        self, path, header=True, index=True,
        sep='\t', float_format='%.1f'
    ):
        self.data.to_csv(
            path, header=header, index=index,
            sep=sep, float_format=float_format)

    def view(self):
        view = TableView()
        model = DataFrameModel(self.data)
        view.setModel(model)
        return view


class RecordTotal(RecordTable):
    export_name = 'summary_total.tsv'

    def export(self, path: Path) -> None:
        self._export_table(path, header=False)

    def view(self):
        view = TableView()
        model = SeriesModel(self.data)
        view.setModel(model)
        view.horizontalHeader().hide()
        return view


class RecordByTaxon(RecordTable):
    export_name = 'summary_per_sample.tsv'

    def export(self, path: Path) -> None:
        self._export_table(path, index=False)


class RecordByGene(RecordTable):
    export_name = 'summary_per_marker.tsv'

    def export(self, path: Path) -> None:
        self._export_table(path)


class RecordByInput(RecordTable):
    export_name = 'summary_per_input_file.tsv'

    def export(self, path: Path) -> None:
        self._export_table(path)


class RecordDisjoint(RecordData):
    export_name = 'disjoint_groups.txt'

    def export(self, path: Path) -> None:
        with open(path, 'w') as file:
            for count, group in enumerate(self.data, start=1):
                print(f'Group {count}', file=file)
                print('-------------', file=file)
                for sample in group:
                    print(sample, file=file)
                print('', file=file)

    def view(self):
        view = TreeView()
        model = DisjointGroupsModel(self.data)
        view.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        view.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        view.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        view.setModel(model)
        view.expandAll()
        view.header().hide()
        return view


class RecordForeign(RecordData):
    export_name = 'no_overlap_pairs.txt'

    def export(self, path: Path) -> None:
        with open(path, 'w') as file:
            for x, y in self.data:
                print(f'{x}\t{y}', file=file)

    def view(self):
        view = TableView()
        model = ForeignPairsModel(self.data)
        view.setModel(model)
        return view


class RecordOutliers(RecordData):
    export_name = 'outliers.txt'

    def export(self, path: Path) -> None:
        with open(path, 'w') as file:
            for gene, data in self.data.items():
                if not data:
                    continue
                print(f'{gene}', file=file)
                print('-------------', file=file)
                if isinstance(data, str):
                    print(data, file=file)
                elif isinstance(data, list):
                    for sample in data:
                        print(sample, file=file)
                print('', file=file)

    def view(self):
        view = TreeView()
        model = OutliersModel(self.data)
        view.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        view.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        view.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        view.setModel(model)
        view.expandAll()
        view.header().hide()
        return view


class RecordPadded(RecordTable):
    export_name = 'padded.tsv'

    def export(self, path: Path) -> None:
        self._export_table(path)


class SummaryReport:
    def __init__(
        self,
        total: Record,
        by_taxon: Record,
        by_gene: Record,
        by_input: Record,
    ):
        self.records = AttrDict()
        self.records.total = total
        self.records.by_taxon = by_taxon
        self.records.by_gene = by_gene
        self.records.by_input = by_input


class Diagnoser:
    def __init__(self, params: DiagnoserParams):
        self.params = params
        self.strings = self._get_strings()
        self.op_general_info = OpGeneralInfo()
        self.op_general_info_per_gene = OpGeneralInfoPerGene()
        self.op_general_info_per_file = OpGeneralInfoPerFile()
        self.op_sequence_bouncer = OpSequenceBouncer(self.params.iqr)
        self.temp = TemporaryDirectory(prefix='concat_diagnose_')
        self.filename = None
        self.timestamp = None
        self.scheme = None
        self.input_files_count = None
        self.disjoint_groups = None
        self.foreign_pairs = None
        self.outlier_count = None
        self._writer_padding = False
        self._writer_frames = False

    def pipe_input_streams(self, input_streams):
        self.input_files_count = len(input_streams)
        if not self.params.report:
            return input_streams
        return [
            stream.pipe(self.op_general_info_per_file)
            for stream in input_streams]

    def pipe_aligned_stream(self, aligned_stream):
        if not self.params.report:
            return aligned_stream
        return aligned_stream.pipe(OpGeneralInfoTagMafftRealigned())

    def _pipe_padding_info(self, stream):
        if not self.params.report:
            return stream
        if self._writer_padding:
            stream = stream.pipe(OpGeneralInfoTagPaddedLength())
        if self._writer_frames:
            stream = stream.pipe(OpGeneralInfoTagPaddedCodonPosition())
        return stream

    def _pipe_general_info(self, stream):
        if not any([
            self.params.report,
            self.params.disjoint,
            self.params.foreign
        ]):
            return stream
        return (
            stream
            .pipe(self.op_general_info)
            .pipe(self.op_general_info_per_gene)
            )

    def _pipe_sequence_bouncer(self, stream):
        if not self.params.outliers:
            return stream
        return stream.pipe(self.op_sequence_bouncer)

    def modify_writer_filters(self, writer):
        self._writer_padding = bool(getattr(writer, 'padding', False))
        self._writer_frames = bool(getattr(writer, 'adjust_frames', False))
        writer.filters = [
            self._pipe_padding_info,
            writer.filter,
            self._pipe_general_info,
            self._pipe_sequence_bouncer,
        ]

    def _get_table_header(self):
        return pd.Series({
            'Output filename': self.filename,
            'Output timestamp': self.timestamp,
            'Output format': self.scheme.text,
            'Number of input files used': self.input_files_count,
        })

    def get_dataframe(self):
        return self.op_general_info.table.dataframe

    def _get_table_total(self):
        header = self._get_table_header()
        table = self.op_general_info.table.total_data()
        table["GC content (%)"] = f'{table["GC content (%)"]:.1f}'
        table["Missing nucleotides (%)"] = f'{table["Missing nucleotides (%)"]:.1f}'
        table["Average number of nucleotides per sample"] = f'{table["Average number of nucleotides per sample"]:.2f}'
        table["Average number of markers per sample"] = f'{table["Average number of markers per sample"]:.2f}'
        table["Average number of samples per marker"] = f'{table["Average number of samples per marker"]:.2f}'
        return pd.concat([header, table])

    def _get_table_by_taxon(self):
        return self.op_general_info.table.by_taxon()

    def _get_table_by_gene(self):
        table = self.op_general_info.table
        return self.op_general_info_per_gene.get_info(table)

    def _get_table_by_input_file(self):
        return self.op_general_info_per_file.get_info()

    def get_summary_report(self) -> Optional[SummaryReport]:
        strings = self.strings['report']
        formatter = self._export_formatter()
        if self.params.report:
            return SummaryReport(
                total = Record(
                    RecordFlag.Void,
                    strings['total']['title'],
                    strings['total']['description'],
                    RecordTotal(self._get_table_total(), formatter)),
                by_taxon = Record(
                    RecordFlag.Void,
                    strings['by_taxon']['title'],
                    strings['by_taxon']['description'],
                    RecordByTaxon(self._get_table_by_taxon(), formatter)),
                by_gene = Record(
                    RecordFlag.Void,
                    strings['by_gene']['title'],
                    strings['by_gene']['description'],
                    RecordByGene(self._get_table_by_gene(), formatter)),
                by_input = Record(
                    RecordFlag.Void,
                    strings['by_input']['title'],
                    strings['by_input']['description'],
                    RecordByInput(self._get_table_by_input_file(), formatter)),
            )
        return None

    def _get_strings(self):
        path = resources.get('itaxotools.concatenator_gui', 'docs/diagnoser.json')
        with open(path) as file:
            return json.load(file)

    def _get_record_disjoint(self) -> Record:
        if not self.params.disjoint:
            return None
        strings = self.strings['disjoint']
        groups = list(list(group) for group in self.op_general_info.table.disjoint_taxon_groups())
        if len(groups) <= 1:
            return Record(
                RecordFlag.Info,
                strings['Info']['title'],
                strings['Info']['description'])
        else:
            formatter = self._export_formatter()
            return Record(
                RecordFlag.Warn,
                strings['Warn']['title'],
                strings['Warn']['description'],
                data=RecordDisjoint(groups, formatter))

    def _get_record_foreign(self) -> Record:
        if not self.params.foreign:
            return None
        strings = self.strings['foreign']
        pairs = list(self.op_general_info.table.unconnected_taxons())
        if not pairs:
            return Record(
                RecordFlag.Info,
                strings['Info']['title'],
                strings['Info']['description'])
        else:
            formatter = self._export_formatter()
            return Record(
                RecordFlag.Warn,
                strings['Warn']['title'],
                strings['Warn']['description'],
                data=RecordForeign(pairs, formatter))

    def _get_record_outliers(self) -> Record:
        if not self.params.outliers:
            return None
        strings = self.strings['outliers']
        outliers = self.op_sequence_bouncer.outliers
        count = sum(
            len(data) for gene, data in outliers.items()
            if isinstance(data, list))
        if not count:
            return Record(
                RecordFlag.Info,
                strings['Info']['title'],
                strings['Info']['description'])
        else:
            formatter = self._export_formatter()
            return Record(
                RecordFlag.Warn,
                strings['Warn']['title'],
                strings['Warn']['description'],
                data=RecordOutliers(outliers, formatter))

    def _get_record_padded(self) -> Record:
        if not self.params.report:
            return None
        strings = self.strings['padded']
        table = self._get_table_by_gene()
        padded = table['Padded to compensate for unequal sequence lengths'] == 'yes'
        if any(padded):
            formatter = self._export_formatter()
            return Record(
                RecordFlag.Warn,
                strings['Warn']['title'],
                strings['Warn']['description'],
                data=RecordPadded(table[padded][[
                    'Padded to compensate for unequal sequence lengths',
                    'Re-aligned by MAFFT',
                    ]], formatter))
        return None

    def get_record_log(self) -> RecordLog:
        log = RecordLog()
        log.append(self._get_record_padded())
        log.append(self._get_record_disjoint())
        log.append(self._get_record_foreign())
        log.append(self._get_record_outliers())
        return log

    def export_dir(self):
        return Path(self.temp.name)

    def _export_record(self, record):
        path = self.export_dir() / record.data.export_name
        record.data.export(path)

    def export_all(self):
        temp = self.export_dir()
        if self.params.report:
            for record in self.get_summary_report().records.values():
                self._export_record(record)
        for record in self.get_record_log():
            if record.data is not None:
                self._export_record(record)

    def _export_formatter(self):
        if self.filename is None:
            return None
        return f'{self.filename.stem}_{{}}'


class SummaryReportLabel(QtWidgets.QLabel):
    clicked = QtCore.Signal(str)

    def __init__(self, attr, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setCursor(QtCore.Qt.PointingHandCursor)
        self.attr = attr

    def enterEvent(self, event):
        font = self.font()
        font.setUnderline(True)
        self.setFont(font)

    def leaveEvent(self, event):
        font = self.font()
        font.setUnderline(False)
        self.setFont(font)

    def mousePressEvent(self, event):
        if (
            event.type() == QtCore.QEvent.MouseButtonPress and
            event.button() == QtCore.Qt.LeftButton
        ):
            self.clicked.emit(self.attr)


class SummaryReportView(QtWidgets.QWidget):
    link_clicked = QtCore.Signal(object)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.report = None

        self.total = SummaryReportLabel('total', 'Total')
        self.per_input = SummaryReportLabel('by_input', 'per input file')
        self.per_sample = SummaryReportLabel('by_taxon', 'per sample')
        self.per_marker = SummaryReportLabel('by_gene', 'per marker')

        self.total.clicked.connect(self._clicked)
        self.per_input.clicked.connect(self._clicked)
        self.per_sample.clicked.connect(self._clicked)
        self.per_marker.clicked.connect(self._clicked)

        layout = QtWidgets.QHBoxLayout()
        layout.addSpacing(10)
        layout.addWidget(QtWidgets.QLabel('\u2b95'))
        layout.addSpacing(6)
        layout.addWidget(QtWidgets.QLabel('Summary reports:'))
        layout.addSpacing(8)
        layout.addWidget(self.total)
        layout.addWidget(QtWidgets.QLabel(', '))
        layout.addWidget(self.per_input)
        layout.addWidget(QtWidgets.QLabel(', '))
        layout.addWidget(self.per_sample)
        layout.addWidget(QtWidgets.QLabel(', '))
        layout.addWidget(self.per_marker)
        layout.addStretch(1)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.setLayout(layout)

    def _clicked(self, attr):
        if self.report:
            self.link_clicked.emit(getattr(self.report.records, attr))

    def setReport(self, report: Optional[SummaryReport]):
        self.report = report
        self.setVisible(report is not None)
