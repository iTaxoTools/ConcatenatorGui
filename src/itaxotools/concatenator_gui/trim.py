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

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from itaxotools.concatenator.library.model import GeneSeries, Operator
from itaxotools.pygblocks import compute_mask, trim_sequence


class OpTrimGblocks(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.genes = set()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        print(f'Trimming {gene.name}:\n')
        if gene.series is None:
            print('N/A')
            print()
            return gene

        print('- Original sequences:\n')
        for sequence in gene.series:
            print(sequence)
        print()

        print('- Mask:\n')
        mask = compute_mask(sequence for sequence in gene.series)
        print(mask)
        print()

        gene.series = gene.series.apply(trim_sequence, mask=mask)
        print('- Trimmed sequences:\n')
        for sequence in gene.series:
            print(sequence)

        print(f'\n{"-"*20}\n')

        self.genes.add(gene.name)
        return gene


class OpTrimClipKit(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.genes = set()

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:

        msa = MultipleSeqAlignment(
            SeqRecord(Seq(seq), id=id)
            for id, seq in gene.series.items())

        print(f'Trimming {gene.name}:\n')
        print(msa)
        print(f'\n{"-"*20}\n')

        self.genes.add(gene.name)
        return gene

