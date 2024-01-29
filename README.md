# ConcatenatorGui

Perform a sequence of transformations on a series of input files:

* Concatenate sequences from different files
* Convert between different file formats: Nexus, Fasta, Tabfile etc.
* Rename or delete markers and create codon subsets.
* Align sequences per marker using MAFFTpy.
* Trim sequences per marker using pyGblocks or ClipKIT.
* Calculate phylogenetic trees per marker using FastTree.

This is a Qt GUI for [Concatenator](https://github.com/iTaxoTools/concatenator).


### Windows and macOS Executables
Download and run the standalone executables without installing Python.</br>
[See the latest release here.](https://github.com/iTaxoTools/ConcatenatorGui/releases/latest)


### Installing from source
Clone and install the latest version (requires Python 3.8.6 or later):
```
git clone https://github.com/iTaxoTools/ConcatenatorGui.git
cd ConcatenatorGui
pip install . -f packages.html
```


## Usage
To launch the GUI, please use:
```
concatenator-gui
```

The program will then guide you through a series of steps.

It will usually be possible to hover the available options or information with the mouse cursor to display additional help in the form of tooltips.

It is also possible to sort the displayed markers by a specific field (name, nucleotides etc.) by clicking on the corresponding column header.


#### 1. Import Input Files

Import one or more sequence files by drag-and-drop or by using the "Import" button. Imported files are briefly checked and some quick statistics of their content will be displayed on the screen. You may expand a file to inspect the contained markers by double-clicking it.


#### 2. Filter Genes

You may rename or delete the imported markers. Gene names must be unique.


#### 3. Align Sequences

You may choose to align your sequences per marker using [MAFFTpy](https://github.com/iTaxoTools/MAFFTpy). First select a strategy, which will be used for all alignments. Then highlight the markers you wish to align and click "Align" to select them. Finally, click "Next" to begin the calculations.


#### 4. Trim Sequences

You may choose to trim your sequences per marker using [pyGblocks](https://github.com/iTaxoTools/pygblocks) or [ClipKIT](https://github.com/JLSteenwyk/ClipKIT). First select one of either toolkits, which will be used for all alignments. A set of parameters will appear depending on the tool selected, which can be customized or kept at the default values. Continue by clicking "Next", then highlight the markers you wish to trim and click "Trim" to select them. Finally, click "Next" to begin the calculations.


#### 5. Character Sets by Codon

You may subset any number of markers by codon. This information is only exported for certain file formats, such as Nexus. It is possible to specify the reading frame of each marker and the names of the resulting character sets. For bulk editing, select one or more markers and click the "Edit" button. You may also double-click a specific field on the marker table to edit it.

In the future, it will be possible to automatically check for invalid reading frames and to determine the reading frames for each marker based on the genetic code type.


#### 6. Export Sequence Data

Start by selecting an output file format, then set the corresponding options. You may hover any option with the mouse cursor for more information. When you are ready, click "Export" and select a save location.

You may additionally calculate phylogenetic trees using [FastTreePy](https://github.com/iTaxoTools/FastTreePy). This option is only available if all markers are of the same length. You may either calculate one tree for the whole alignment or a tree for each marker. The trees will be saved at the same location as the sequence data.

Optionally, you may perform data validation tests on the dataset and produce summary reports per input file, per marker, per sample and for the whole dataset. This includes the option to exhaustively detect outlier sequences per marker using [SequenceBouncer](https://github.com/iTaxoTools/SequenceBouncer).


### Supported Formats

The following sequence file formats are available:

* Interleaved Nexus
* Concatenated Fasta, Phylip or Ali
* Multi-file Fasta, Phylip or Ali
* PartitionFinder and IQTree
* Tab-separated vector files

File formats with multiple files may be directly imported/exported in the form of either directories or zip archives.

Tabfiles must follow a specific format. The first column must be the "species", optionally followed by descriptor columns such as "voucher" or "locality". Each marker column must have a name prefixed with "sequences_" (eg. "sequence_gene01"). See the provided example file [here](examples/test_tabfile.tab).


### Packaging

It is recommended to use PyInstaller from within a virtual environment:
```
pip install ".[dev]" -f packages.html
pyinstaller scripts/concatenator.spec
```


### Citing

Concatenator was developed by V. Kharchev and S. Patmanidis in the framework of the iTaxoTools project:

Vences, M., A. Miralles, S. Brouillet, J. Ducasse, A. Fedosov, V. Kharchev, I. Kostadinov, S. Kumari, S. Patmanidis, M.D. Scherz, N. Puillandre, S.S. Renner (2021):
iTaxoTools 0.1: Kickstarting a specimen-based software toolkit for taxonomists. - Megataxa 6: 77-92.

Concatenator integrates Mafft (by K. Katoh and collaborators), FastTree (by M. N. Price) and SequenceBouncer (by C. D. Dunn).
