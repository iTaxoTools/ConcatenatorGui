
from dataclasses import dataclass, field
from pathlib import Path
from typing import Set, Optional


@dataclass(frozen=True)
class Sample:
    name: str


@dataclass
class Charset:
    name: str
    nucleotides: int = 0
    missing: int = 0
    uniform: bool = False
    samples: Set[Sample] = field(default_factory=set, repr=False)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)


@dataclass
class File:
    path: Path
    format: Optional[str] = None
    nucleotides: int = 0
    missing: int = 0
    uniform: bool = False
    samples: Set[Sample] = field(default_factory=set, repr=False)
    charsets: Set[Charset] = field(default_factory=set, repr=False)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.path == other.path

    def __hash__(self):
        return hash(self.path)


class Concatenation:
    def __init__(self):
        self.files = set()
        self.charsets = set()
        self.samples = set()
