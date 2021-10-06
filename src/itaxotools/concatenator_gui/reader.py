
from typing import Callable, Dict
from pathlib import Path

import pandas as pd

from itaxotools.concatenator.library.file_types import FileType
from itaxotools.concatenator.library.detect_file_type import autodetect
from itaxotools.concatenator.library.nexus import read as nexus_read


CallableReader = Callable[[Path], pd.DataFrame]

type_reader: Dict[FileType, CallableReader] = dict()


def reader(type: FileType) -> CallableReader:
    def decorator(func: CallableReader) -> CallableReader:
        type_reader[type] = func
        return func
    return decorator


@reader(FileType.TabFile)
def readTabFile(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep='\t')
    return data


@reader(FileType.NexusFile)
def readNexusFile(path: Path) -> pd.DataFrame:
    with path.open() as file:
        data = nexus_read(file)
    return data


class ReaderNotFound(Exception):
    def __init__(self, type: FileType):
        self.type = type
        super().__init__(f'No reader for FileType: {str(type)}')


def dataframe_from_path(path: Path) -> pd.DataFrame:
    type = autodetect(path)
    if not type in type_reader:
        raise ReaderNotFound(type)
    data = type_reader[type](path)
    return data
