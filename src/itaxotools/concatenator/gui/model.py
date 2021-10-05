
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Set, Optional


@dataclass(frozen=True)
class Sample:
    name: str


@dataclass
class Charset:
    name: str  # key
    nucleotides: int = 0
    missing: int = 0
    uniform: bool = False
    samples: Set[Sample] = field(default_factory=set, repr=False)


@dataclass
class File:
    path: Path  # key
    format: Optional[str] = None
    nucleotides: int = 0
    missing: int = 0
    uniform: bool = False
    samples: Set[Sample] = field(default_factory=set, repr=False)
    charsets: Dict[str, Charset] = field(default_factory=dict, repr=False)

    @property
    def name(self):
        return self.path.name

class Concatenation:
    def __init__(self):
        self.files: Dict[Path, File] = dict()
        self.charsets: Dict[str, Charset] = dict()
        self.samples: Set[Sample] = set()
