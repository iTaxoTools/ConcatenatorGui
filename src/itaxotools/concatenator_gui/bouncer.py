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

"""SequenceBouncer operator"""


from typing import Dict, List, Union
from logging import Handler, Formatter, INFO
from time import sleep

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from itaxotools.concatenator.library.model import Operator
from itaxotools.concatenator.library.utils import Field

from itaxotools.sequence_bouncer import SequenceBouncer, InputSequence
from itaxotools.sequence_bouncer.SequenceBouncer import logger, version


class PrintHandler(Handler):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setLevel(INFO)

    def emit(self, message):
        print(self.format(message))


class NoLogSequenceBouncer(SequenceBouncer):

    def initialize_logger(self):
        "Also adds file handler if required"
        v = self.vars
        p = self.params

        self.logger = logger.getChild(str(id(self)))
        self.logger.propagate = False
        self.logger.addHandler(PrintHandler())

        self.logger.info('\nSequenceBouncer: A method to remove outlier entries from a multiple sequence alignment\n')
        self.logger.info('Cory Dunn')
        self.logger.info('University of Helsinki')
        self.logger.info('cory.dunn@helsinki.fi')
        self.logger.info('Version: ' + version)
        self.logger.info('Please cite DOI: 10.1101/2020.11.24.395459')
        self.logger.info('___\n')


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
        bouncer = NoLogSequenceBouncer(
            input, IQR_coefficient=self.iqr, gap_percent_cut=100, write_none=True)
        try:
            verdict = bouncer()
            self.outliers[gene.name] = [
                sample for sample, accepted in verdict.items() if not accepted]
        except Exception as exception:
            self.outliers[gene.name] = str(exception)
        return gene
