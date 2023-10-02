**Table of content**
- [Description](#description)
- [Package installation](#package-installation)
- [Usage](#usage)
    - [1.1. Retrieve fasta files from a general database](#11-retrieve-fasta-files-from-a-general-database)
    - [1.2. Create the initial database](#12-create-the-initial-database)
    - [1.3. Filter the sequences](#13-filter-the-sequences)
    - [1.4. Create SA taxons](#14-create-sa-taxons)
  - [2. Merge databases](#2-merge-databases)
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
ASVmaker is designed to be used by modules : ![Drag Racing](res/usage_steps.png)


### 1.1. Retrieve fasta files from a general database
The first step involves downloading a FASTA file for a specific genus of interest from a general database : Silva, Unite, RNAcentral, ENA, NCBI or DDBJ. This file contains the genomic data necessary for subsequent analysis. 


```shell
# The genus of interest
GENUS=$1
# The name of the database. Pick among : ncbi,ebi,ddbj,rnaCentral,silva,unite
DB=$2
# Fasta file to proceed
SEQ=$3
```


### 1.2. Create the initial database
Next, ASVmaker enables the creation of a genus-specific database by using the downloaded FASTA file. Each sequence lineage is verified by accession number through the European Nucleotide Archive API if possible and through the NCBI Entrez API if the ENA one doesn't match. Users must specify primers to be used during the simulation of the amplification process, allowing for precise targeting of the desired genomic regions and ASV creation.

```shell
#  INIT DATABASE
python3 -m asvmaker \
    -inf database/${DB}/${GENUS}_${DB}_info_create.txt \
    create \
    -i ${SEQ} \
    -db ${DB} \
    -fp fw_primer.fasta \
    -rp rv_primer.fasta \
    -fmt 5 \
    -rmt 5 \
    -o database/${DB}/${GENUS}_${DB}_create.json
```
### 1.3. Filter the sequences
To enhance the quality and specificity of the analysis, ASVmaker provides the functionality to filter out redundant amplicons and exclude unwanted taxonomy. Redundant amplicons are ASVs sharing the same taxonomy. Unwanted taxonomy or species that are not of interest (e.g. : “sp.” or “aff.”) can also be filtered out, ensuring a more focused analysis on the target genus. 

```shell
#  FILTER DATABASE
python3 -m asvmaker \
    -inf database/${DB}/${GENUS}_${DB}_info_filter.txt \
    filter \
    -i database/${DB}/${GENUS}_${DB}_create.json \
    -g1 ${GENUS} \
    -o database/${DB}/${GENUS}_${DB}_filter.json

#  EXPORT RESULTS
python3 -m asvmaker \
    -inf database/${DB}/${GENUS}_${DB}_info_exp1.txt \
    export \
    -i database/${DB}/${GENUS}_${DB}_filter.json \
    -sao database/${DB}/${GENUS}_${DB}_sa.txt 1
```

### 1.4. Create SA taxons

ASVmaker creates Shared Amplicon (SA) groups, which involve clustering identical ASVs which have different taxonomies. This grouping allows for a comprehensive understanding of the taxonomic diversity within the selected ASV, providing valuable precisions into the composition and dynamics of microbial communities.

```shell
#  EDIT DATABASE
python3 -m asvmaker \
    -inf database/${DB}/${GENUS}_${DB}_info_edit.txt \
    edit \
    -i database/${DB}/${GENUS}_${DB}_filter.json \
    -grp database/${DB}/${GENUS}_${DB}_sa_ext.txt \
    -o database/${DB}/${GENUS}_${DB}_edit.json

#  EXPORT RESULTS
python3 -m asvmaker \
    -inf database/${DB}/${GENUS}_${DB}_info_exp2.txt \
    export \
    -i database/${DB}/${GENUS}_${DB}_edit.json \
    -aop database/${DB}/${GENUS}_${DB}_asv.fasta
```

## 2. Merge databases

Moreover, ASVmaker offers the option to merge ASV specific databases from different general databases, providing flexibility to combine data from various sources. This merging process allows for a more comprehensive dataset, enabling comparative analysis and broader insights into the studied genus.

```shell
#  MERGE DATABASES
python3 -m asvmaker \
    -inf database/all/${GENUS}_info_merge.txt \
    merge \
    -i database/rnaCentral/${GENUS}_rnaCentral_edit.json \
    -sa1 database/rnaCentral/${GENUS}_rnaCentral_sa_ext.txt \
    -i2 database/unite/${GENUS}_unite_edit.json \
    -sa2 database/unite/${GENUS}_unite_sa_ext.txt \
    -o database/all/${GENUS}_merge.json

#  EXPORT RESULTS
python3 -m asvmaker \
    -inf database/all/${GENUS}_info_exp3.txt \
    export \
    -i database/all/${GENUS}_merge.json \
    -aop database/all/${GENUS}_asv.fasta
```


# References
Article in revision. 