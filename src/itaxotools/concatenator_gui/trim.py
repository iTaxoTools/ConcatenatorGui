# -----------------------------------------------------------------------------
# ConcatenatorQt - GUI for Concatenator
# Copyright (C) 2021-2023  Patmanidis Stefanos
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
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
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

        mode = TrimmingMode.smart_gap

        alignment = MultipleSeqAlignment(
            SeqRecord(Seq(seq), id=id)
            for id, seq in gene.series.items())

        sequence_type = get_seq_type(alignment)
        gap_characters = get_gap_chars(sequence_type)
        gaps = smart_gap_threshold_determination(alignment, gap_characters)

        msa = create_msa(alignment, gap_characters)
        msa.trim(mode, gap_threshold=gaps)

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
