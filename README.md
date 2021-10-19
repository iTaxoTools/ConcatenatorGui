# ConcatenatorGui

This is a Qt GUI for Concatenator: <https://github.com/iTaxoTools/concatenator>

Perform a sequence of transformations on a series of input files:

* Convert between different file formats: Nexus, Fasta, Tabfile etc.
* Rename or delete character sets and create codon subsets.
* Align sequences per character set using MAFFT.

## Quick start

```
$ pip install . -f packages.html
$ concatenator-gui examples/test_tabfile.tab
```

## Packaging

```
$ pip install ".[dev]" -f packages.html
$ pyinstaller scripts/concatenator.spec
```
