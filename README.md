**Table of content**
- [Description](#description)
- [Package installation](#package-installation)
- [Usage](#usage)
  - [Create database](#create-database)
  - [Merge databases](#merge-databases)
- [References](#references)

# Description

This package allow to create ASV specific reference databases from unfiltered FASTA files. FASTA file supported:
- [ncbi: 'nt'](https://www.ncbi.nlm.nih.gov/nucleotide/)
- [ebi: 'ena'](https://www.ebi.ac.uk/ena/browser/home)
- [ddbj: 'arsa'](http://ddbj.nig.ac.jp/arsa/)
- [RNAcentral](https://rnacentral.org/)
- [UNITE](https://unite.ut.ee/repository.php)
- [Silva](https://www.arb-silva.de/)

It is mainly intended to be used by command line. However, some modules can be used inside python scripts. Most of the database creation is automatically done, only a few commands are left to the user.

CAUTION: this tool is not a pipeline but a set of modules.

# Package installation

The package is available on  [pypi](https://pypi.org/). It can be installed by the following command: 

```shell
pip install asvmaker
```

# Usage

## Create database

1) Retrieve fasta files from a general database
2) Create the initial database
3) Filter the sequences
4) Create SA taxons


## Merge databases


# References