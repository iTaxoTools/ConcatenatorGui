
from __future__ import annotations

from typing import Dict, List, Set, Hashable, Iterable, Optional
from dataclasses import dataclass, field
from pathlib import Path


class DataSet:
    def __init__(self):
        self.data: Dict[Hashable, int] = dict()
        self.groups: List[DataGroup] = list()
        self.index_lengths: Dict[int, int] = dict()
        self.index_count: int = 0

    def index_next(self) -> int:
        self.index_count += 1
        return self.index_count

    def update(self, index: int, items: Iterable[Hashable]) -> int:
        if index not in self.index_lengths:
            self.data.update({item: index for item in items})
            self.index_lengths[index] = len(items)
            return index
        if len(items) == self.index_lengths[index]:
            return index
        else:
            new_index = self.index_next()
            for item in items:
                self.data[item] = new_index
            groups = [g for g in self.groups if index in g.indices]
            for group in groups:
                group.indices.add(new_index)
            return new_index

    def remove(self, group: DataGroup):
        ...


class DataGroup:
    def __init__(self, dataset: DataSet, items: Iterable[Hashable] = []):
        self.dataset: DataSet = dataset
        self.indices: Set[int] = set()
        self.length: int = 0
        self.index: int = dataset.index_next()
        self.update(items)
        dataset.groups.append(self)

    def update(self, items: Iterable[Hashable]):
        queue: Dict[int, List[Hashable]] = dict()

        for item in items:
            if item in self.dataset.data:
                if not self.dataset.data[item] in self.indices:
                    self.length += 1
                index = self.dataset.data[item]
            else:
                self.length += 1
                index = self.index
            if not index in queue:
                queue[index] = list()
            queue[index].append(item)

        for index in queue:
            real_index = self.dataset.update(index, queue[index])
            self.indices.add(real_index)

    def __len__(self):
        return self.length

    @property
    def data(self):
        ...


@dataclass(frozen=True)
class Sample:
    name: str


@dataclass
class Charset:
    name: str  # key
    nucleotides: int = 0
    missing: int = 0
    uniform: bool = False
    samples: Optional[DataGroup] = field(default=None, repr=False)


@dataclass
class File:
    path: Path  # key
    format: Optional[str] = None
    nucleotides: int = 0
    missing: int = 0
    uniform: bool = False
    samples: Optional[DataGroup] = field(default=None, repr=False)
    charsets: Dict[str, Charset] = field(default_factory=dict, repr=False)

    @property
    def name(self):
        return self.path.name

class Concatenation:
    def __init__(self):
        self.files: Dict[Path, File] = dict()
        self.charsets: Dict[str, Charset] = dict()
        self.samples: DataSet = set()
