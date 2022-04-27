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

"""SequenceBouncer operator"""


from typing import Dict, List, Union
from time import sleep

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from itaxotools.concatenator.library.model import Operator
from itaxotools.concatenator.library.utils import Field

from itaxotools.sequence_bouncer import SequenceBouncer, InputSequence


class OpSequenceBouncer(Operator):
    iqr: float = Field('iqr', value=1.0)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.outliers: Dict[str, Union[str, List[str]]] = dict()

    def input_from_gene(self, gene):
        iterator = iter(
            SeqRecord(Seq(sequence), id=str(name), name=str(name))
            for name, sequence in gene.series.items()
        )
        return InputSequence(gene.name, iterator)

    def call(self, gene):
        input = self.input_from_gene(gene)
        bouncer = SequenceBouncer(
            input, IQR_coefficient=self.iqr, gap_percent_cut=100, write_none=True)
        try:
            verdict = bouncer()
            self.outliers[gene.name] = [
                sample for sample, accepted in verdict.items() if not accepted]
        except Exception as exception:
            self.outliers[gene.name] = str(exception)
        print(f'BOUNCED {gene.name}', self.outliers[gene.name])
        return gene
