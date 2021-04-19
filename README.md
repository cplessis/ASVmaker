# Table of contents


- [Description](#description)
- [Project Structure](#project-structure)
- [Package installation](#package-installation)
  + [From Python Index (Pypi)](#from-python-index--pypi-)
  + [Conda environment](#conda-environment)
  + [Python library directory](#python-library-directory)
- [Command lines](#command-lines)
  + [Calling the package](#calling-the-package)
  + [Required parameters](#required-parameters)
  + [Optionals parameters](#optionals-parameters)
- [Module for import usage](#module-for-import-usage)



# Description

This package is used for creating reference databases from unfiltered FASTA files. The FASTA files should come from a general purpose database such as [ncbi: 'nt'](https://www.ncbi.nlm.nih.gov/nucleotide/), [ebi: 'ena'](https://www.ebi.ac.uk/ena/browser/home), or [ddbj: 'arsa'](http://ddbj.nig.ac.jp/arsa/). And the format of the input FASTA files must be as following (complete description +  single or multiple lines sequence) :

```
>JN176092.1 Fusarium mexicanum strain MXJAL-19 translation elongation factor 1 alpha (EF1) gene, partial cds
GTCGACTCTGGCAAGTCGACCACTGTGAGTACAACCCTCGACGATGAGCTTATCTGCCATCGTCATCCCG
ATCACCATCGATATTGCTCTCTGGAAGTTCGAGACTCCTCGCTACTATGTCACCGTCATTGGTATGTTGT
CGCTCATGCCTCGTTCTCCCTTTATTCGTACTAACATATCACTCAGACGCTCCC
```

It is mainly intended to be used with command lines with python. However, some modules can be used for other python projects by import. Most of the database creation is automatically done, only a few commands are left to the user. After the creation of the files, a manual verification should be done. A usage example is available [here](https://github.com/cplessis/Fusarium-EF1A-Q2_RefDb) in which we create a reference database for the  Fusarium EF1 alpha gene. 



# Project Structure

``` shell
.
├── LICENCE.txt
├── README.md
├── pyproject.toml
├── setup.cfg
├── src
│   └── q2_mkrefdb
│       ├── __init__.py
│       ├── __main__.py
│       ├── fastafilter.py
│       ├── pcr.py
│       └── utils.py
└── tests
```

* [LICENCE](./LICENCE.txt) :: MIT
* [pyproject](./pyproject.toml) :: Use for Pypi
* [setup](./setup.cfg) :: Use for Pypi
* [src/q2_mkrefdb](./src/q2_mkrefdb) :: Package
  * [init](./src/q2_mkrefdb/__init__.py) :: Initiate Package when imported
  * [main](./src/q2_mkrefdb/__main__.py) :: Used for command lines
  * [fastafilter](./src/q2_mkrefdb/fastafilter.py) :: Main algorithm
  * [pcr](./src/q2_mkrefdb/pcr.py) :: Used in fastafliter
  * [utils](./src/q2_mkrefdb/utils.py) :: Used in pcr and fastafilter
* [tests](./tests) :: Validation tests



# Package installation

### From Python Index (Pypi)

The package is available on  [pypi](https://pypi.org/). This installation method is the easiest one. 

```shell
pip install q2_mkrefdb
```

OR

```shell
python -m pip install q2_mkrefdb
```

### Conda environment

Installing a package directly with Pypi can cause some version conflicts. Also, the use of this package is generally linked to the use of [Qiime2](https://qiime2.org/). You will therefore find a Conda environment [here](https://github.com/cplessis/Fusarium-EF1A-Q2_RefDb) which will contain all the necessary packages as well as the version of Qiime2 which corresponds to the analyses we carried out in the tutorial [Fusarium-EF1A-Q2_RefDb](https://github.com/cplessis/Fusarium-EF1A-Q2_RefDb).

### Python library directory

If you would like to install the package offline, you can use this slightly less "conventional" method. Start by downloading the "q2_mkrefdb" directory from [src](./src). Once you have the directory (package), place it in your python library with the other packages you already use with python. Usually a large part of the python packages are in the "site-packages" directory. You can find the location of this directory by running the following command with `python` (where `anymodule` is the name of the module of your choice):

```python
import anymodule
print(anymodule.__file__)
```

Once the package is placed in the library directory you can call it like any other python package.

# Commands lines

### Calling the package

``` shell
python -m q2_mkrefb [args]
```

### Required parameters

```shell
-i, --sequences_input SEQUENCES_INPUT
# File path to the raw sequences FASTA file to treat
"The file format must be FASTA. The sequences can be written on one or multiple line. The sequences description must start with '>' and must be followed by the accession number. Such as : '>KJ679381|KJ679381.1 Fusarium falciforme ...' OR '>ENA|KJ679394|KJ679394.1 Fusarium keratoplasticum ...' OR '>JN176092.1 Fusarium mexicanum ...'. These three examples are respectively from DDBJ, EBI, NCBI databases.
```



 ```shell
-db, --source_database {ncbi,ebi,ddbj}
# Name of the database from which the FASTA file was obtained.
 ```



```shell
-fp, --forward_primer FORWARD_PRIMER
# File path to the forward primer FASTA file.
```

 

```shell
-rp, --reverse_primer REVERSE_PRIMER
# File path to the reverse primer FASTA file.
```

 

### Optionals parameters

```shell
-fmt, --fw_mismatch_tol FW_MISMATCH_TOL
# Number of mismatch which are accepted for the forward primer annealing.
```

  

```shell
-rmt, --rv_mismatch_tol RV_MISMATCH_TOL
# Number of mismatch which are accepted for the reverse primer annealing.
```

 

```shell
-inf, --infos_file INFOS_FILE
# Informations file in which all the files treatments are recorded.
```



```shell
-dit, --displ_inf_terminal
# Display information file in terminal if arg is specified.
```

 

```shell
-g1, --genus1 GENUS1
# Name of the first genus to keep. Ex: Fusarium
```

 

```shell
-g2, --genus2 GENUS2
# Name of the second genus to keep. Ex: Gibberella
```

 

```shell
-u, --unverified
# Keep the unverified genus if arg is specified.
```

 

```shell
-tbc, --taxonomy_by_complex
# Do not group species with same complex lineage if arg is specified.
```

 

```shell
-f, --filtering_type {taxon,global}
# Filtering by taxon or globally
```

 

```shell
-rsv, --redund_seq_variants
# Keep the redundant sequences variants if arg is specified.
```

 

```shell
-ra, --redundant_amplicon
# Keep the redundant amplicons if arg is specified.
```

 

```shell
-nasv, --seq_without_ampl
# Keep the sequences without amplicons if arg is specified.
```

---



```shell
-svo, --seq_variants_output SEQ_VARIANTS_OUTPUT
# File name of the sequences variants FASTA output.
```

   

```shell
-ao, --amplicons_output AMPLICONS_OUTPUT
# File name of the amplicons FASTA output.
```

 

```shell
-to TAXONOMY_OUTPUT, --taxonomy_output TAXONOMY_OUTPUT
# File name of the taxonomy TXT output.
```

 

```shell
-tlo, --taxon_list_output TAXON_LIST_OUTPUT
# File name of the taxon list TXT output.
```

 

```shell
-sao, --shared_ampl_output SHARED_AMPL_OUTPUT
# File name of the shared amplicons TXT output.
```

 

```shell
-cdo, --complex_dict_output COMPLEX_DICT_OUTPUT
# File name of the complex dictionnary TXT output.
```

  

```shell
-mto, --modified_tax_output MODIFIED_TAX_OUTPUT
# File name of the modified taxon list TXT output.
```

 

```shell
-ado, --access_dict_output ACCESS_DICT_OUTPUT
# File name of the access dictionnary JSON output.
```

  

# Module for import usage

More informations will be available soon...



