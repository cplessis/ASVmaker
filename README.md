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
    + [Filtration](#filtration)
    + [Export](#export)
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
├── LICENSE.txt
├── README.md
├── pyproject.toml
├── setup.cfg
├── src
│   └── q2_mkrefdb
│       ├── __init__.py
│       ├── __main__.py
│       ├── database.py
│       ├── edit.py
│       ├── filter.py
│       ├── lineage.py
│       ├── merge.py
│       ├── pcr.py
│       ├── stats.py
│       └── utils.py
└── tests
```

* [LICENSE](./LICENSE.txt) :: MIT
* [pyproject](./pyproject.toml) :: Use for Pypi
* [setup](./setup.cfg) :: Use for Pypi
* [src/q2_mkrefdb](./src/q2_mkrefdb) :: Package
  * [init](./src/q2_mkrefdb/__init__.py) :: Initiate Package when imported
  * [main](./src/q2_mkrefdb/__main__.py) :: Use for command lines
  * [database](./src/q2_mkrefdb/database.py) :: Main algorithm
  * [edit](./src/q2_mkrefdb/edit.py) :: 
  * [filter](./src/q2_mkrefdb/filter.py) :: 
  * [lineage](./src/q2_mkrefdb/lineage.py) ::
  * [merge](./src/q2_mkrefdb/merge.py) ::
  * [pcr](./src/q2_mkrefdb/pcr.py) :: Used in fastafliter
  * [stats](./src/q2_mkrefdb/stats.py) ::
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

Once the package has been placed in the library directory you can call it like any other python package.

# Command lines

## CREATE

``` shell
python -m q2_mkrefb create [args]
```

### Required

---

```shell
-i, --sequences_input SEQUENCES_INPUT
# File path to the raw sequences FASTA file to treat.
```

The file format must be FASTA. The sequences can be written on one or multiple lines. The sequences description must start with '>' and must be followed by the accession number. Such as : `>KJ679381|KJ679381.1 Fusarium falciforme ...`OR `>ENA|KJ679394|KJ679394.1 Fusarium keratoplasticum ...` OR  `>JN176092.1 Fusarium mexicanum ...`. These three examples are respectively from DDBJ, EBI, NCBI databases. We highly recommand to use the description format of one of these database in order to avoid any error.

---

 ```shell
-db, --source_database {ncbi,ebi,ddbj}
# Name of the database from which the FASTA file was obtained.
 ```

Each of the databases (ncbi, ebi, or ddjb) have its own format of description for the fasta files. **NCBI**: `>JN176092.1 Fusarium mexicanum ...`; **EBI**: `>ENA|KJ679394|KJ679394.1 Fusarium keratoplasticum ...`; **DDBJ**: `>KJ679381|KJ679381.1 Fusarium falciforme ...`. If you are not using one of these databases, choose the one wich have the same sequences description format.

---

```shell
-fp, --forward_primer FORWARD_PRIMER
-rp, --reverse_primer REVERSE_PRIMER
# File path to the forward primer FASTA file.
```

The primer FASTA file must have only one sequence. With a description and a sequence.  You need one FASTA file per primer.

----

```shell
-fmt, --fw_mismatch_tol FW_MISMATCH_TOL
-rmt, --rv_mismatch_tol RV_MISMATCH_TOL
# Number of mismatch which are accepted for the forward ('-fmt') and reverse ('-rmt') primer annealing.
```

 The mismatch tolerance is the number of mismatch you accept on the primer annealing site. The default value is 3 per primer. You must consider IUPAC nucleotide code as a mismatch. Our algorithm does NOT take into account the degenerated nucloetides yet. A primer can only anneal to a site where the primer 3'-end is perflectly matching to the binding sequence.

---

```shell
-tp, --trim_primers
# This option keep the primers on amplicons if specified. By default the primers are trimmed.
```

 By default, the primers are removed from the amplicon. If this option is specified, they will be conserved. 

### Optional

---

```shell
-inf, --infos_file INFOS_FILE
# Informations file in which all the files treatments are recorded.
```

The INFOS_FILE save all the database information (number of sequences, amplicons, etc.) at every step of the analysis. As soon as the data are modified by a filtration, the actual state of the data is written on the file.

___

```shell
-dit, --displ_inf_terminal
# Display information file in terminal if arg is specified.
```

 If you specify this parameter, all the database states will be printed in your shell during the analysis. This command can be specified or not independantbly from the `--infos_file` arg. 

---

```shell
-mto, --modified_tax_output MODIFIED_TAX_OUTPUT
# File name of the modified taxon list TXT output.
```

Export all the modified taxon in file with as first column the old name of taxon and as second column the new name of the taxon after lineage verification. 

___

___

## FILTER

```shell
python -m q2_mkrefb filter [args]
```

### Required

---

```shell
-i, --database_json DATABASE_JSON
# Path to the database JSON file to treat.
```

The inputed DATABASE_JSON file must have been created with the CREATE option. This file contain all the informations about the database which will be filtered. It is the output of the -o, --output_database_json parameter. 

---

```shell
-o, --output_database_json OUTPUT_DATABASE_JSON
# File name of the access dictionnary JSON output.
```

Export the whole dictionnary of the database from the programm. This dictionnary is the last state of the database before exporting the files. In this file you will find all the information the database. The dictionnary is made as : {accession_number: {taxon: 'the genus_specie', sequence: 'full sequence', lineage: 'lineage', amplicon: 'amplicon sequence', description: 'the sequence description from the original fasta file' }}. The AccessDictionnary can be usefull to analyse statistics about the data.

### Optional

---

```shell
-g1, --genus1 GENUS1
-g2, --genus2 GENUS2
# Name of the first and second genus to keep. Ex: -g1 Fusarium -g2 Gibberella
```

You can specify up to 2 genus to keep in your dataset. All the species which have a different Genus will be removed from the dataset. Be aware to write the genus with the right spelling and the word must start with an Uppercase. 

---

```shell
-cc, --custom_complex
# File path to the custom complex CSV file.
```

The custom complex arg should be used to create customed complex/group. You need a custom complex CSV file presente as : ***column 1*** = 'new complex name' ***column 2*** = 'taxon name' . All the sequences with the 'taxon name' will be grouped as 'new complex name'. Then all the filtration steps will occur.

---

```shell
-u, --unverified
# Keep the unverified taxon if arg is specified.
```

Some species are not well identified and can be classed as "Fusarium sp." for example. If you specify the `-u` arg these unverified species will be conserved. By default, the unverified species are removed. 

---

```shell
-tbc, --taxonomy_by_complex
# Do not group species with same complex lineage if arg is specified.
```

 By default the species are grouped by complex. A complex is a group of closely related species. This information about species is taken from the lineage.  If you specify the `-tbc` arg the species will NOT be group by complex.

---

```shell
-rsv, --redund_seq_variants
# Keep the redundant sequences variants if arg is specified.
```

 By default the redundant sequences variants are removed from the dataset. If this arg is specified, they are conserved. 

---

```shell
-ra, --redundant_amplicon
# Keep the redundant amplicons if arg is specified.
```

  By default the redundant amplicons are removed from the dataset. If this arg is specified, they are conserved. 

---

```shell
-nasv, --seq_without_ampl
# Keep the sequences without amplicons if arg is specified.
```

By default the sequences which can not amplify the region of interest are removed from the dataset. If this arg is specified they are conserved.

---

```shell
-inf, --infos_file INFOS_FILE
# Informations file in which all the files treatments are recorded.
```

The INFOS_FILE save all the database information (number of sequences, amplicons, etc.) at every step of the analysis. As soon as the data are modified by a filtration, the actual state of the data is written on the file.

___

```shell
-dit, --displ_inf_terminal
# Display information file in terminal if arg is specified.
```

 If you specify this parameter, all the database states will be printed in your shell during the analysis. This command can be specified or not independantbly from the `--infos_file` arg. 

___

___

## EDIT

```shell
python -m q2_mkrefb edit [args]
```

### Required

---

```shell
-i, --database_json DATABASE_JSON
# Path to the database JSON file to treat.
```

The inputed DATABASE_JSON file must have been created with the CREATE option. This file contain all the informations about the database which will be filtered. It is the output of the -o, --output_database_json parameter. 

---

```shell
-o, --output_database_json OUTPUT_DATABASE_JSON
# File name of the access dictionnary JSON output.
```

Export the whole dictionnary of the database from the programm. This dictionnary is the last state of the database before exporting the files. In this file you will find all the information the database. The dictionnary is made as : {accession_number: {taxon: 'the genus_specie', sequence: 'full sequence', lineage: 'lineage', amplicon: 'amplicon sequence', description: 'the sequence description from the original fasta file' }}. The AccessDictionnary can be usefull to analyse statistics about the data.

### Optional

---

```shell
-rm, --remove REMOVE
# Remove all the sequences with IDs on id_list_csv.
```

The CSV file must be as : ***column 1*** = 'seq_id' 

---

```shell
-mv, --rename RENAME
# Rename all the sequences with IDs on id_list_csv.
```

The CSV file must be as : ***column 1*** = 'seq_id'  ***column 2*** = 'new_taxon_name'

---

```shell
-grp, --group GROUP
# Path to the database JSON file to treat.
```

Group all the sequences with IDs on shared_ext_csv on a commune taxon name. The taxon name will become the SA_id of the group. Only one amplicon will represent the group after this command is runned. The shared_ext_csv FILE must be the one generated with EXPORT shared ampl.

---

```shell
-inf, --infos_file INFOS_FILE
# Informations file in which all the files treatments are recorded.
```

The INFOS_FILE save all the database information (number of sequences, amplicons, etc.) at every step of the analysis. As soon as the data are modified by a filtration, the actual state of the data is written on the file.

___

```shell
-dit, --displ_inf_terminal
# Display information file in terminal if arg is specified.
```

 If you specify this parameter, all the database states will be printed in your shell during the analysis. This command can be specified or not independantbly from the `--infos_file` arg. 

___

___

## EXPORT

```shell
python -m q2_mkrefb export [args]
```

### Required

---

```shell
-i, --database_json DATABASE_JSON
# Path to the database JSON file to treat.
```

The inputed DATABASE_JSON file must have been created with the CREATE option. This file contain all the informations about the database which will be filtered. It is the output of the -o, --output_database_json parameter. 

### Optional

____

```shell
-svo, --seq_variants_output SEQ_VARIANTS_OUTPUT
# File name of the sequences variants FASTA output.
```

Export a file with all the sequences variants as FASTA format with the related accession number as description. These are the full length sequences after all filtrations. One sequence is related to one amplicon.

---

```shell
-ao, --amplicons_output AMPLICONS_OUTPUT
# File name of the amplicons FASTA output.
```

Export a file with all the amplicons as FASTA format with the related accession number as description. One amplicon is related to one full sequence. 

---

```shell
-aop, --ampl_output_phylo AMPL_OUTPUT_PHYLO
# File name of the amplicon variants FASTA output in NCBI format (for phylogeny).
```

Export a file with all the amplicons as FASTA format with the related accession number and taxon as description. This file should be used for phylogeny ([https://ngphylogeny.fr/](https://ngphylogeny.fr/)).

---

```shell
-to, --taxonomy_output TAXONOMY_OUTPUT
# File name of the taxonomy TXT output.
```

Export a file with 2 columns. The first contain the accession number (sequence id) and the second contain the related lineage of this id. 

---

```shell
-tlo, --taxon_list_output TAXON_LIST_OUTPUT
# File name of the taxon list TXT output.
```

Export a file with the list of all the taxon after filtration. One line, one taxon. 

---

```shell
-sao, --shared_ampl_output SHARED_AMPL_OUTPUT
# File name of the shared amplicons TXT output.
```

Export a file with the sequences which have exactly the same amplicon. One line is related to one amplicon. For example, if three sequence id are written on one line, that means their amplicons are all exactly the same. 

---

```shell
-cdo, --complex_dict_output COMPLEX_DICT_OUTPUT
# File name of the complex dictionnary TXT output.
```

Export  a file with the list of all the complex names as first column and then a dictionnary as second column with all the accession number (related to ASV) and species of this complex.

---

```shell
-inf, --infos_file INFOS_FILE
# Informations file in which all the files treatments are recorded.
```

The INFOS_FILE save all the database information (number of sequences, amplicons, etc.) at every step of the analysis. As soon as the data are modified by a filtration, the actual state of the data is written on the file.

___

```shell
-dit, --displ_inf_terminal
# Display information file in terminal if arg is specified.
```

 If you specify this parameter, all the database states will be printed in your shell during the analysis. This command can be specified or not independantbly from the `--infos_file` arg. 

# Module for import usage

More informations will be available soon...



