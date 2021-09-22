# ConcatenatorGui

This is a Qt GUI for Concatenator: <https://github.com/iTaxoTools/concatenator>

Perform a sequence of transformations on a series of input files:

* Convert between different file formats
* Rename or delete character sets
* Align sequences per character set
* Subset character sets by codon

_This GUI is a development prototype and is not yet connected to the Concatenator core:_

* _Random contents will be generated for each file and action._
* _No changes will be written to disk._

## Quick start

```
$ pip install .
$ concatenator-gui
```

## Packaging

It is advised to do this inside a virtual environment using a tool such as pipenv:

```
$ pipenv shell
$ pip install -e .[dev]
$ pyinstaller scripts/concatenator.spec
```
