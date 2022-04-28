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


from tempfile import TemporaryDirectory
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from itaxotools.common.utility import AttrDict
from itaxotools.concatenator.library.operators import (
    OpChainGenes, OpTranslateGenes, OpApplyToGene, OpTagSet,
    OpUpdateMetadata, OpFilterGenes, OpGeneralInfo,
    OpGeneralInfoPerFile, OpGeneralInfoPerGene,
    OpGeneralInfoTagMafftRealigned,
    OpGeneralInfoTagPaddedLength,
    OpGeneralInfoTagPaddedCodonPosition)

from .bouncer import OpSequenceBouncer
from .records import Record, RecordFlag, RecordLog


@dataclass
class DiagnoserParams:
    report: bool = True
    disjoint: bool = True
    foreign: bool = True
    outliers: bool = False
    iqr: float = 20.0


class SummaryReport:
    def __init__(
        self,
        total: pd.DataFrame,
        by_taxon: pd.DataFrame,
        by_gene: pd.DataFrame,
        by_input: pd.DataFrame,
    ):
        self.records = AttrDict()
        self.records.total = Record(RecordFlag.Info, 'Total', data=total)
        self.records.by_taxon = Record(RecordFlag.Info, 'Per Taxon', data=by_taxon)
        self.records.by_gene = Record(RecordFlag.Info, 'Per Gene', data=by_gene)
        self.records.by_input = Record(RecordFlag.Info, 'Per Input File', data=by_input)


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

    def get_table_header(self):
        return pd.Series({
            'Name of output file': self.filename,
            'Timestamp of output file': self.timestamp,
            'Output file format': self.scheme.text,
            'Number of input files used': self.input_files_count,
        })

    def get_dataframe(self):
        return self.op_general_info.table.dataframe

    def _get_table_total(self):
        header = self.get_table_header()
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

    def get_summary_report(self) -> SummaryReport:
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
            return Record(RecordFlag.Warn, 'Disjoint WARN', data=groups)

    def _get_record_foreign(self) -> Record:
        if not self.params.foreign:
            return None
        pairs = list(self.op_general_info.table.unconnected_taxons())
        if not pairs:
            return Record(RecordFlag.Info, 'Foreign OK')
        else:
            return Record(RecordFlag.Warn, 'Foreign WARN', data=pairs)

    def _get_record_outliers(self) -> Record:
        if not self.params.outliers:
            return None
        outliers = self.op_sequence_bouncer.outliers
        count = sum(len(data) for gene, data in outliers.items() if isinstance(data, list))
        if not count:
            return Record(RecordFlag.Info, 'Outliers OK')
        else:
            return Record(RecordFlag.Warn, 'Outliers WARN', data=outliers)

    def get_record_log(self) -> RecordLog:
        log = RecordLog()
        log.append(self._get_record_disjoint())
        log.append(self._get_record_foreign())
        log.append(self._get_record_outliers())
        return log

    def export_dir(self):
        return Path(self.temp.name)

    def export_table(
        self, table, name, header=True, index=True,
        sep='\t', float_format='%.1f'
    ):
        table.to_csv(
            self.export_dir() / name, header=header, index=index,
            sep=sep, float_format=float_format)

    def export_disjoint(self, name):
        self.disjoint_groups = 0
        with open(self.export_dir() / name, 'w') as file:
            for group in self.op_general_info.table.disjoint_taxon_groups():
                self.disjoint_groups += 1
                print(f'Group {self.disjoint_groups}', file=file)
                print('-------------', file=file)
                for sample in group:
                    print(sample, file=file)
                print('', file=file)

    def export_foreign(self, name):
        self.foreign_pairs = 0
        with open(self.export_dir() / name, 'w') as file:
            for x, y in self.op_general_info.table.unconnected_taxons():
                self.foreign_pairs += 1
                print(f'{x}\t{y}', file=file)
            if self.foreign_pairs == 0:
                print('No foreign pairs detected.', file=file)

    def export_outliers(self, name):
        self.outlier_count = 0
        with open(self.export_dir() / name, 'w') as file:
            for gene, data in self.op_sequence_bouncer.outliers.items():
                if not data:
                    continue
                print(f'{gene}', file=file)
                print('-------------', file=file)
                if isinstance(data, str):
                    print(data, file=file)
                elif isinstance(data, list):
                    self.outlier_count += len(data)
                    for sample in data:
                        print(sample, file=file)
                print('', file=file)
            if self.outlier_count == 0:
                print('No outliers detected.', file=file)

    def export_all(self):
        temp = self.export_dir()
        if self.params.report:
            self.export_table(self._get_table_total(), 'total_data.tsv', header=False)
            self.export_table(self._get_table_by_input_file(), 'by_input_file.tsv', index=False)
            self.export_table(self._get_table_by_taxon(), 'by_taxon.tsv')
            self.export_table(self._get_table_by_gene(), 'by_gene.tsv')
        if self.params.disjoint:
            self.export_disjoint('disjoint_groups.txt')
        if self.params.foreign:
            self.export_foreign('foreign_pairs.txt')
        if self.params.outliers:
            self.export_outliers('outliers.txt')
