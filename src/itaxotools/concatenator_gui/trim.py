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

"""Trimming operators"""

from typing import Optional
from enum import Enum

import pandas as pd

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from itaxotools.concatenator.library.model import GeneSeries, Operator
from itaxotools.pygblocks import Options, compute_mask, trim_sequence

from clipkit.helpers import get_seq_type, get_gap_chars, create_msa
from clipkit.smart_gap_helper import smart_gap_threshold_determination
from clipkit.modes import TrimmingMode
from clipkit.msa import MSA


class ClipKitTrimmingMode(Enum):
    smart_gap = "smart-gap", "dynamic determination of gaps threshold"
    gappy = "gappy", "trim all sites that are above a threshold of gappyness"
    kpic = "kpic", "keep only parismony informative and constant sites"
    kpic_smart_gap = "kpic-smart-gap", "a combination of kpic- and smart-gap-based trimming"
    kpic_gappy = "kpic-gappy", "a combination of kpic- and gappy-based trimming"
    kpi = "kpi", "keep only parsimony informative sites"
    kpi_smart_gap = "kpi-smart-gap", "a combination of kpi- and smart-gap-based trimming"
    kpi_gappy = "kpi-gappy", "a combination of kpi- and gappy-based trimming"

    def __init__(self, mode: str, description: str):
        self.mode = mode
        self.description = description


class OpTrimGblocks(Operator):
    def __init__(self, options: Options, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.options = options
        self.genes = set()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        options = Options(**self.options.as_dict())

        print(f'Trimming {gene.name}:\n')
        if gene.series is None:
            print('N/A')
            print()
            return gene

        id_length = max(len(index) for index in gene.series.index) + 4

        for id, sequence in gene.series.items():
            print(id.ljust(id_length), sequence)
        print()

        mask = compute_mask((sequence for sequence in gene.series), options)
        print('MASK'.ljust(id_length), mask)
        print()

        gene.series = gene.series.apply(trim_sequence, mask=mask)
        for id, sequence in gene.series.items():
            print(id.ljust(id_length), sequence)

        print(f'\n{"-"*20}\n')

        self.genes.add(gene.name)
        return gene


class OpTrimClipKit(Operator):
    def __init__(self, options: dict, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mode = getattr(TrimmingMode, options.get('mode'))
        self.gap_characters = list(options.get('gap_characters', None))
        self.gaps = options.get('gaps', None)
        self.genes = set()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        print(f'Trimming {gene.name}:\n')
        if gene.series is None:
            print('N/A')
            print()
            return gene

        id_length = max(len(index) for index in gene.series.index) + 4

        for id, sequence in gene.series.items():
            print(id.ljust(id_length), sequence)
        print()

        alignment = MultipleSeqAlignment(
            SeqRecord(Seq(seq), id=id)
            for id, seq in gene.series.items())

        if not self.gap_characters:
            sequence_type = get_seq_type(alignment)
            self.gap_characters = get_gap_chars(sequence_type)

        if self.mode in {
            TrimmingMode.smart_gap,
            TrimmingMode.kpi_smart_gap,
            TrimmingMode.kpic_smart_gap,
        }:
            self.gaps = smart_gap_threshold_determination(alignment, self.gap_characters)

        msa = create_msa(alignment, self.gap_characters)
        msa.trim(self.mode, gap_threshold=self.gaps)

        for index, record in zip(gene.series.index, msa.to_bio_msa()):
            gene.series.at[index] = str(record.seq)

        mask = self.get_mask(msa)
        print('MASK'.ljust(id_length), mask)
        print()

        for id, sequence in gene.series.items():
            print(id.ljust(id_length), sequence)

        print(f'\n{"-"*20}\n')

        self.genes.add(gene.name)
        return gene

    def get_mask(self, msa: MSA) -> str:
        length = msa._original_length
        trimmed = msa._site_positions_to_trim
        return ''.join('.' if x in trimmed else '#' for x in range(length))
