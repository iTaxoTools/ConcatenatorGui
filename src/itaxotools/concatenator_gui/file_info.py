
from pathlib import Path

from itaxotools.concatenator.library.detect_file_type import autodetect
from itaxotools.concatenator.library.file_types import FileType

from itaxotools.concatenator_gui.file_iterator import iterator_from_path

from . import model


type_short = {
    FileType.TabFile: 'Tabfile',
    FileType.NexusFile: 'Nexus',
    FileType.FastaFile: 'Fasta',
    FileType.PhylipFile: 'Phylip',
    FileType.ConcatTabFile: 'Tabfile',
    FileType.ConcatFasta: 'Fasta',
    FileType.ConcatPhylip: 'Phylip',
    FileType.MultiFastaOutput: 'Fasta (zip)',
    FileType.MultiPhylipOutput: 'Phylip (zip)',
    FileType.MultiAliOutput: 'Ali (zip)',
    FileType.MultiFastaInput: 'Fasta (zip)',
    FileType.MultiPhylipInput: 'Phylip (zip)',
    FileType.MultiAliInput: 'Ali (zip)',
    FileType.PartitionFinderOutput: 'PartitionFinder',
    FileType.CodonTab: 'CodonTab',
}


def file_info_from_path(
    path: Path,
    samples: model.DataSet
) -> model.File:

    type = autodetect(path)
    data = iterator_from_path(path)
    file = model.File(path)
    all_characters = 0
    all_characters_missing = 0
    all_uniform = []

    for series in data:
        seq = series.name
        lengths = series.str.len()
        missing = series.str.count('-')
        len_test = len(series[0])
        uniform = all(series.str.len() == len_test)
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

    file.format = type_short[type]
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
