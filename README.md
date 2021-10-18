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

It is advised to do this inside a virtual environment using a tool such as pipenv:

```
$ pipenv shell
$ pip install -e ".[dev]"
$ pyinstaller scripts/concatenator.spec
```
