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

"""Get simple file info from path"""


from typing import Callable
from pathlib import Path

import re

from itaxotools.concatenator import (
    FileType, FileFormat, autodetect, read_from_path)
from itaxotools.concatenator.library.utils import has_uniform_length

from . import model


formats_short = {
    FileType.File: {
        FileFormat.Tab: 'Tabfile',
        FileFormat.Nexus: 'Nexus',
        FileFormat.Ali: 'Ali',
        FileFormat.Fasta: 'Fasta',
        FileFormat.Phylip: 'Phylip',
    },
    FileType.Directory: {
        FileFormat.Ali: 'Multi-Ali',
        FileFormat.Fasta: 'Multi-Fasta',
        FileFormat.Phylip: 'Multi-Phylip',
    },
    FileType.ZipArchive: {
        FileFormat.Ali: 'Multi-Ali',
        FileFormat.Fasta: 'Multi-Fasta',
        FileFormat.Phylip: 'Multi-Phylip',
    },
}


def file_info_from_path(
    path: Path,
    samples: model.DataSet,
    checker: Callable[[None], None] = None,
) -> model.File:

    type, format = autodetect(path)
    stream = read_from_path(path)
    file = model.File(path)
    all_characters = 0
    all_characters_missing = 0
    all_uniform = []

    for gene in stream:
        series = gene.series
        seq = series.name
        lengths = series.str.len()
        missing = series.str.count(f'[{re.escape(gene.missing + gene.gap)}]')
        uniform = has_uniform_length(series)
        mask = (lengths - missing != 0)
        samples_new = series.index[mask]

        charset = model.Charset(seq)
        charset.characters = sum(lengths)
        charset.characters_missing = sum(missing)
        charset.uniform = 'Yes' if uniform else 'No'
        charset.samples = model.DataGroup(samples)
        charset.samples.update(samples_new)

        file.charsets[seq] = charset
        all_characters += charset.characters
        all_characters_missing += charset.characters_missing
        all_uniform += [uniform]
        if checker is not None:
            checker()

    file.format = formats_short[type][format]
    file.characters = all_characters
    file.characters_missing = all_characters_missing
    if all(all_uniform):
        file.uniform = 'Yes'
    elif all(not x for x in all_uniform):
        file.uniform = 'No'
    else:
        file.uniform = 'Mixed'
    file.samples = model.DataGroup(samples)
    file.samples.merge(cs.samples for cs in file.charsets.values())
    return file
