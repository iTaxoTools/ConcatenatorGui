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

"""Data diagnoser"""


from typing import Optional
from tempfile import TemporaryDirectory
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from PySide6.QtCore import Qt, QSize, QRect, QEvent, QAbstractItemModel, QAbstractListModel, Signal
from PySide6.QtWidgets import QWidget, QListView, QStyledItemDelegate, QStyle, QSizePolicy, QHBoxLayout, QLabel
from PySide6.QtGui import QFontMetrics

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
    outliers: bool = True
    iqr: float = 20.0


class RecordTable(RecordData):
    def _export_table(
        self, path, header=True, index=True,
        sep='\t', float_format='%.1f'
    ):
        self.data.to_csv(
            path, header=header, index=index,
            sep=sep, float_format=float_format)


class RecordTotal(RecordTable):
    export_name = 'total_data.tsv'

    def export(self, path: Path) -> None:
        self._export_table(path, header=False)


class RecordByTaxon(RecordTable):
    export_name = 'by_taxon.tsv'

    def export(self, path: Path) -> None:
        self._export_table(path, index=False)


class RecordByGene(RecordTable):
    export_name = 'by_gene.tsv'

    def export(self, path: Path) -> None:
        self._export_table(path)


class RecordByInput(RecordTable):
    export_name = 'by_input_file.tsv'

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


class RecordForeign(RecordData):
    export_name = 'foreign_pairs.txt'

    def export(self, path: Path) -> None:
        with open(path, 'w') as file:
            for x, y in self.data:
                print(f'{x}\t{y}', file=file)


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


class SummaryReport:
    def __init__(
        self,
        total: pd.DataFrame,
        by_taxon: pd.DataFrame,
        by_gene: pd.DataFrame,
        by_input: pd.DataFrame,
    ):
        self.records = AttrDict()
        self.records.total = Record(RecordFlag.Info, 'Total', data=RecordTotal(total))
        self.records.by_taxon = Record(RecordFlag.Info, 'Per Taxon', data=RecordByTaxon(by_taxon))
        self.records.by_gene = Record(RecordFlag.Info, 'Per Gene', data=RecordByGene(by_gene))
        self.records.by_input = Record(RecordFlag.Info, 'Per Input File', data=RecordByInput(by_input))


class Diagnoser:
    def __init__(self, params: DiagnoserParams):
        self.params = params
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
            'Name of output file': self.filename,
            'Timestamp of output file': self.timestamp,
            'Output file format': self.scheme.text,
            'Number of input files used': self.input_files_count,
        })

    def get_dataframe(self):
        return self.op_general_info.table.dataframe

    def _get_table_total(self):
        header = self._get_table_header()
        table = self.op_general_info.table.total_data()
        table["GC content of sequences"] = f'{table["GC content of sequences"]:.1f}'
        table["% missing data"] = f'{table["% missing data"]:.1f}'
        table["Average number of nucleotides per taxon"] = f'{table["Average number of nucleotides per taxon"]:.2f}'
        table["Average number of markers per taxon"] = f'{table["Average number of markers per taxon"]:.2f}'
        table["Average number of taxa per marker"] = f'{table["Average number of taxa per marker"]:.2f}'
        return pd.concat([header, table])

    def _get_table_by_taxon(self):
        return self.op_general_info.table.by_taxon()

    def _get_table_by_gene(self):
        table = self.op_general_info.table
        return self.op_general_info_per_gene.get_info(table)

    def _get_table_by_input_file(self):
        return self.op_general_info_per_file.get_info()

    def get_summary_report(self) -> Optional[SummaryReport]:
        if self.params.report:
            return SummaryReport(
                total = self._get_table_total(),
                by_taxon = self._get_table_by_taxon(),
                by_gene = self._get_table_by_gene(),
                by_input = self._get_table_by_input_file(),
            )
        return None

    def _get_record_disjoint(self) -> Record:
        if not self.params.disjoint:
            return None
        # Todo: save in a tree formation
        groups = list(self.op_general_info.table.disjoint_taxon_groups())
        if len(groups) <= 1:
            return Record(RecordFlag.Info, 'Disjoint OK')
        else:
            return Record(RecordFlag.Warn, 'Disjoint WARN', data=RecordDisjoint(groups))

    def _get_record_foreign(self) -> Record:
        if not self.params.foreign:
            return None
        pairs = list(self.op_general_info.table.unconnected_taxons())
        if not pairs:
            return Record(RecordFlag.Info, 'Foreign OK')
        else:
            return Record(RecordFlag.Warn, 'Foreign WARN', data=RecordForeign(pairs))

    def _get_record_outliers(self) -> Record:
        if not self.params.outliers:
            return None
        outliers = self.op_sequence_bouncer.outliers
        count = sum(len(data) for gene, data in outliers.items() if isinstance(data, list))
        if not count:
            return Record(RecordFlag.Info, 'Outliers OK')
        else:
            return Record(RecordFlag.Warn, 'Outliers WARN', data=RecordOutliers(outliers))

    def get_record_log(self) -> RecordLog:
        log = RecordLog()
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
            print(str(record))
            if record.data is not None:
                self._export_record(record)


class SummaryReportLabel(QLabel):
    clicked = Signal(str)

    def __init__(self, attr, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setCursor(Qt.PointingHandCursor)
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
            event.type() == QEvent.MouseButtonPress and
            event.button() == Qt.LeftButton
        ):
            self.clicked.emit(self.attr)


class SummaryReportView(QWidget):
    clicked = Signal(object)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.report = None

        self.total = SummaryReportLabel('total', 'Total')
        self.per_input = SummaryReportLabel('by_input', 'per input')
        self.per_sample = SummaryReportLabel('by_taxon', 'per sample')
        self.per_marker = SummaryReportLabel('by_gene', 'per marker')

        self.total.clicked.connect(self._clicked)
        self.per_input.clicked.connect(self._clicked)
        self.per_sample.clicked.connect(self._clicked)
        self.per_marker.clicked.connect(self._clicked)

        layout = QHBoxLayout()
        layout.addSpacing(10)
        layout.addWidget(QLabel('\u2b95'))
        layout.addSpacing(6)
        layout.addWidget(QLabel('Summary reports:'))
        layout.addSpacing(8)
        layout.addWidget(self.total)
        layout.addWidget(QLabel(', '))
        layout.addWidget(self.per_input)
        layout.addWidget(QLabel(', '))
        layout.addWidget(self.per_sample)
        layout.addWidget(QLabel(', '))
        layout.addWidget(self.per_marker)
        layout.addStretch(1)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.setLayout(layout)

    def _clicked(self, attr):
        if self.report:
            self.clicked.emit(getattr(self.report.records, attr))

    def setReport(self, report: SummaryReport):
        self.report = report
