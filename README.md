# ConcatenatorGui

Perform a sequence of transformations on a series of input files:

* Concatenate sequences from different files
* Convert between different file formats: Nexus, Fasta, Tabfile etc.
* Rename or delete genes and create codon subsets.
* Align sequences per gene using MAFFT.
* Calculate phylogenetic trees per gene using FastTree.

This is a Qt GUI for [Concatenator](https://github.com/iTaxoTools/concatenator).


### Executables
Download and run the standalone executables without installing Python.</br>
[See the latest release here.](https://github.com/iTaxoTools/ConcatenatorGui/releases/latest)


### Installing
Clone and install the latest version (requires Python 3.8 or later):
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

It is also possible to sort the displayed genes by a specific field (name, nucleotides etc.) by clicking on the corresponding column header.


#### 1. Import Input Files

Import one or more sequence files by drag-and-drop or by using the "Import" button. Imported files are briefly checked and some quick statistics of their content will be displayed on the screen. You may expand a file to inspect the contained genes by double-clicking it.


#### 2. Filter Genes

You may rename or delete the imported genes. Gene names must be unique.


#### 3. Align Sequences

You may choose to align your sequences per gene using [MAFFTpy](https://github.com/iTaxoTools/MAFFTpy). First select a strategy, which will be used for all gene alignments. Then select the genes you wish to align and click "Align" to mark them. Finally, click "Next" to begin the calculations.


#### 4. Character Sets by Codon

You may subset any number of genes by codon. This information is only exported for certain file formats, such as Nexus. It is possible to specify the reading frame of each gene and the names of the resulting character sets. For bulk editing, select one or more genes and click the "Edit" button. You may also double-click a specific field on the gene table to edit it.

In the future, it will be possible to automatically check for invalid reading frames and to determine the reading frames for each gene based on the genetic code type.


#### 5. Export Sequence Data

Start by selecting an output file format, then set the corresponding options. You may hover any option with the mouse cursor for more information. When you are ready, click "Export" and select a save location.

You may additionally calculate phylogenetic trees using [FastTreePy](https://github.com/iTaxoTools/FastTreePy). This option is only available if all genes are of the same length. You may either calculate one tree for the whole alignment or a tree for each gene. The trees will be saved at the same location as the sequence data.


### Supported Formats

The following sequence file formats are available:

* Interleaved Nexus
* Concatenated Fasta, Phylip or Ali
* Multi-file Fasta, Phylip or Ali
* PartitionFinder and IQTree
* Tab-separated vector files

File formats with multiple files may be directly imported/exported in the form of either directories or zip archives.

Tabfiles must follow a specific format. The first column must be the "species", optionally followed by descriptor columns such as "voucher" or "locality". Each gene column must have a name prefixed with "sequences_" (eg. "sequence_gene01"). See the provided example file [here](examples/test_tabfile.tab).


### Packaging

It is recommended to use PyInstaller from within a virtual environment:
```
pip install ".[dev]" -f packages.html
pyinstaller scripts/concatenator.spec
```
