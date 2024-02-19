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

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, NamedTuple
from collections import defaultdict

import re

from itaxotools.concatenator.library.model import Operator, GeneSeries
from itaxotools.concatenator.library.utils import Field


@dataclass
class ExclusionParams:
    by_site: bool = False
    by_marker: bool = False
    maximum_missing_sites: float = 90.00
    maximum_missing_markers: float = 90.00


class MissingInfo(NamedTuple):
    length: int
    missing: int


class OpCountMissing(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data: dict[str, MissingInfo] = defaultdict(list)

    def call(self, gene: GeneSeries) -> Optional[GeneSeries]:
        lengths = gene.series.str.len()
        missing = gene.series.str.count(f'[{re.escape(gene.missing + gene.gap)}]')
        for (id1, length), (id2, missing) in zip(lengths.items(), missing.items()):
            assert id1 == id2
            info = MissingInfo(length, missing)
            self.data[id1].append(info)
        return gene

    def _should_exclude_id_by_site(self, id: str, maximum_missing_sites: float) -> bool:
        total_length = sum(x.length for x in self.data[id])
        total_missing = sum(x.missing for x in self.data[id])
        if not total_length:
            return True
        if total_missing / total_length >= maximum_missing_sites:
            return True
        return False

    def exclude_by_site(self, maximum_missing_sites: float) -> set[str]:
        return {id for id in self.data if self._should_exclude_id_by_site(id, maximum_missing_sites)}

    def _should_exclude_id_by_marker(self, id: str, total_genes: int, maximum_missing_markers: float) -> bool:
        genes = sum(1 for info in self.data[id] if info.missing < info.length)
        if not total_genes:
            return True
        if genes / total_genes <= maximum_missing_markers:
            return True
        return False

    def exclude_by_marker(self, maximum_missing_markers: float) -> set[str]:
        total_genes = max(len(infos) for infos in self.data.values())
        return {id for id in self.data if self._should_exclude_id_by_marker(id, total_genes, maximum_missing_markers)}
