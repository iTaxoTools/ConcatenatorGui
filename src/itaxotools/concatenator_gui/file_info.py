
from pathlib import Path

from itaxotools.concatenator import (
    FileType, FileFormat, autodetect, read_from_path)
from itaxotools.concatenator.library.operators import OpDropEmpty
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
        FileFormat.Ali: 'MultiAli',
        FileFormat.Fasta: 'MultiFasta',
        FileFormat.Phylip: 'MultiPhylip',
    },
    FileType.ZipArchive: {
        FileFormat.Ali: 'MultiAli (zip)',
        FileFormat.Fasta: 'MultiFasta (zip)',
        FileFormat.Phylip: 'MultiPhylip (zip)',
    },
}


def file_info_from_path(
    path: Path,
    samples: model.DataSet
) -> model.File:

    type, format = autodetect(path)
    stream = read_from_path(path)
    file = model.File(path)
    all_characters = 0
    all_characters_missing = 0
    all_uniform = []

    for series in stream:
        series = OpDropEmpty('?-')(series)
        seq = series.name
        lengths = series.str.len()
        missing = series.str.count('[?-]')
        uniform = has_uniform_length(series)
        mask = (lengths - missing != 0)
        species = series.index[mask]

        charset = model.Charset(seq)
        charset.characters = sum(lengths)
        charset.characters_missing = sum(missing)
        charset.uniform = 'Yes' if uniform else 'No'
        charset.samples = model.DataGroup(samples)
        charset.samples.update(species)

        file.charsets[seq] = charset
        all_characters += charset.characters
        all_characters_missing += charset.characters_missing
        all_uniform += [uniform]

    file.format = formats_short[type][format]
    file.characters = all_characters
    file.characters_missing = all_characters_missing
    if all(all_uniform):
        file.uniform = 'Yes'
    elif all([not x for x in all_uniform]):
        file.uniform = 'No'
    else:
        file.uniform = 'Mixed'
    file.samples = model.DataGroup(samples)
    file.samples.merge([cs.samples for cs in file.charsets.values()])
    return file
