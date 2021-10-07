
from typing import Callable, Dict, Iterator
from pathlib import Path

from zipfile import ZipFile
from zipp import Path as ZipPath  # BUGFIX: backport from Python 3.9.1

import pandas as pd

from itaxotools.concatenator.library.file_types import FileType
from itaxotools.concatenator.library.detect_file_type import autodetect

from .file_reader import (
    type_readers, dataframe_from_path,
    readAliSeries, readFastaSeries, readPhylipSeries
    )


CallableIterator = Callable[[Path], Iterator[pd.Series]]

type_iterators: Dict[FileType, CallableIterator] = dict()


def type_iterator(type: FileType) -> CallableIterator:
    def decorator(func: CallableIterator) -> CallableIterator:
        type_iterators[type] = func
        return func
    return decorator


def _iterateZipFile(
    path: Path,
    func: Callable[[Path], pd.Series]
) -> Iterator[pd.Series]:
    archive = ZipFile(path)
    for part in archive.namelist():
        part_path = ZipPath(archive, part)
        yield func(part_path)


@type_iterator(FileType.MultiAliInput)
def iterateMultiAli(path: Path) -> Iterator[pd.Series]:
    for series in _iterateZipFile(path, readAliSeries):
        yield series


@type_iterator(FileType.MultiFastaInput)
def iterateMultiFasta(path: Path) -> Iterator[pd.Series]:
    for series in _iterateZipFile(path, readFastaSeries):
        yield series


@type_iterator(FileType.MultiPhylipInput)
def iterateMultiPhylip(path: Path) -> Iterator[pd.Series]:
    for series in _iterateZipFile(path, readPhylipSeries):
        yield series


class IteratorNotFound(Exception):
    def __init__(self, type: FileType):
        self.type = type
        super().__init__(f'No iterator for FileType: {str(type)}')


def iterator_from_path(path: Path) -> Iterator[pd.Series]:
    """Species as index, sequences as name"""
    type = autodetect(path)
    if type in type_iterators:
        for series in type_iterators[type](path):
            yield series
        return
    if type in type_readers:
        data = dataframe_from_path(path)
        for column in data:
            series = data[column]
            series.name = column
            yield series
        return
    raise IteratorNotFound(type)
