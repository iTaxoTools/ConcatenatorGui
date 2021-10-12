
from __future__ import annotations

from typing import Dict, List, Set, Hashable, Iterable, Optional
from dataclasses import dataclass, field
from pathlib import Path


class DataSet:
    """
    Defines a set of items on which subsets can be defined as DataGroups.
    Efficient when the subsets overlap.
    """
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
        self.groups.remove(group)

    def __len__(self):
        all_indices = set()
        for group in self.groups:
            all_indices.update(group.indices)
        all_data = [x for x in self.data if self.data[x] in all_indices]
        return len(all_data)


class DataGroup:
    def __init__(self, dataset: DataSet, items: Iterable[Hashable] = []):
        self.dataset: DataSet = dataset
        self.indices: Set[int] = set()
        self.length: int = 0
        self.index: int = dataset.index_next()
        self.update(items)
        dataset.groups.append(self)

    def merge(self, groups: Iterable[DataGroup]):
        for group in groups:
            self.indices.update(group.indices)
        self.length = len(self.data)

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
            if index not in queue:
                queue[index] = list()
            queue[index].append(item)

        for index in queue:
            real_index = self.dataset.update(index, queue[index])
            self.indices.add(real_index)

    def __len__(self):
        return self.length

    @property
    def data(self) -> Set[Hashable]:
        data = self.dataset.data
        return [x for x in data if data[x] in self.indices]


@dataclass(frozen=True)
class Sample:
    name: str


@dataclass
class Charset:
    name: str  # key
    characters: int = 0
    characters_missing: int = 0
    uniform: bool = False
    samples: Optional[DataGroup] = field(default=None, repr=False)

    @property
    def nucleotides(self):
        return self.characters - self.characters_missing

    @property
    def missing(self):
        return self.characters_missing / self.characters


@dataclass
class File:
    path: Path  # key
    format: Optional[str] = None
    characters: int = 0
    characters_missing: int = 0
    uniform: bool = False
    samples: Optional[DataGroup] = field(default=None, repr=False)
    charsets: Dict[str, Charset] = field(default_factory=dict, repr=False)

    @property
    def name(self):
        return self.path.name

    @property
    def nucleotides(self):
        return self.characters - self.characters_missing

    @property
    def missing(self):
        return self.characters_missing / self.characters


class Concatenation:
    def __init__(self):
        self.files: Dict[Path, File] = dict()
        self.charsets: Dict[str, Charset] = dict()
        self.samples = DataSet()

    def remove_file(self, file):
        # memory assigned for samples is not cleared
        del self.files[file.path]