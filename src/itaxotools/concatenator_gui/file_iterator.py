
from typing import Callable, Dict, Iterator
from pathlib import Path

import pandas as pd

from itaxotools.concatenator.library.file_types import FileType
from itaxotools.concatenator.library.detect_file_type import autodetect

from .file_reader import type_reader, dataframe_from_path

CallableIterator = Callable[[Path], Iterator[pd.Series]]

type_iterator: Dict[FileType, CallableIterator] = dict()


def iterator(type: FileType) -> CallableIterator:
    def decorator(func: CallableIterator) -> CallableIterator:
        type_iterator[type] = func
        return func
    return decorator


class IteratorNotFound(Exception):
    def __init__(self, type: FileType):
        self.type = type
        super().__init__(f'No iterator for FileType: {str(type)}')


def iterator_from_path(path: Path) -> Iterator[pd.Series]:
    """Species as index, sequences as name"""
    type = autodetect(path)
    if type in type_iterator:
        raise NotImplementedError
    if type in type_reader:
        data = dataframe_from_path(path)
        for column in data:
            series = data[column]
            series.name = column
            yield series
        return
    raise IteratorNotFound(type)
