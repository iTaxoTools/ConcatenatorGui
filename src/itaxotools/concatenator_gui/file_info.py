
from pathlib import Path
from random import randint
from lorem_text import lorem

from itaxotools.concatenator.library.detect_file_type import autodetect
from itaxotools.concatenator.library.file_types import FileType

from itaxotools.concatenator_gui.reader import dataframe_from_path

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
    data = dataframe_from_path(path)
    print(data)
    file = model.File(path)
    file.format = type_short[type]
    file.missing = randint(0, 9999) / 10000
    file.uniform = ['Yes', 'No'][randint(0, 1)]
    for i in range(1, randint(5, 20)):
        name = lorem.words(randint(2, 6)).replace(' ', '_')
        charset = model.Charset(name)
        charset.nucleotides = randint(200, 3000000)
        charset.uniform = ['Yes', 'No'][randint(0, 1)]
        charset.missing = randint(0, 9999) / 10000
        charset.samples = model.DataGroup(samples)
        foobar = range(randint(0, 100), randint(50, 150))
        charset.samples.update(foobar)
        file.charsets[name] = charset
    file.samples = model.DataGroup(samples)
    file.samples.merge([cs.samples for cs in file.charsets.values()])
    return file
