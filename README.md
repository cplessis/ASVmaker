# Description

Ce package permet de créer des bases de données de références à partir de fichiers FASTA non triés issues d'une base de données généraliste telle que [ncbi : 'nt'](https://www.ncbi.nlm.nih.gov/nucleotide/), [ebi : 'ena'](https://www.ebi.ac.uk/ena/browser/home), ou [ddbj : 'arsa](http://ddbj.nig.ac.jp/arsa/). Il est principalement destiné à être utilisé en ligne de commande avec python. Néanmoins certains modules peuvent être utilisé pour d'autres projet python par importation.

La majeur partie du travaille de création de la base de donnée se fait automatiquement, seules quelques commandes sont à la charge de l'utilisateur. Suite à la création des fichiers une vérification manuelle peut être utile. Un exemple d'utilisation est disponible [ici](https://github.com/cplessis/Fusarium-EF1A-Q2_RefDb) dans lequel nous créons une base de données de référence pour le gène EF1 alpha chez Fusarium. 

# Structure du Projet

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



# Installation du package

### Depuis Python Index (Pypi)

Le package est disponible sur [pypi](https://pypi.org/). L'utilisation de cette méthode est la plus classique est la plus simple. Il est peut être directement installé sur votre machine de la manière suivante.

```shell
pip install q2_mkrefdb
```

ou

```shell
python -m pip install q2_mkrefdb
```

### Avec un environnement Conda

l'installation direct d'un package avec Pypi peut poser certains problèmes de versions. Aussi, l'utilisation de ce package est généralement lié à l'utilisation de QIIME2. Vous pourrez donc trouver un environnement Conda [ici]() qui contiendra tous les packages nécessaires ainsi que la version de QIIM2 qui corespond aux analyses que nous avons mené dans le tutorial [Fusarium-EF1A-Q2_RefDb](https://github.com/cplessis/Fusarium-EF1A-Q2_RefDb).

### Dépôt direct du package en librairie

Si vous souhaitez installer le package or connexion, vous pouvez utiliser cette méthode un peu moins "conventionnelle". Commencez par télécharger préalablement le répertoire "q2_mkrefdb" présent dans [src](./src). Une fois le répertoire (package) récupéré, placer dans votre librairie python avec les autres package que vous utilisez déjà en python. Généralement une grande partie des packages python se trouvent dans le répertoire "site-packages". Vous pouvez trouvez la position de ce répertoire en lançant la commande suivante avec python:

```python
import anymodule
print(anymodule.__file__)
```

Une fois la package placé dans le répertoire des librairies vous pourrez l'appeler comme tous les autres packages pythons.

# Lignes de commandes

### Appel du package

``` shell
python -m q2_mkrefb [args]
```

### Arguments obligatoires

```shell
-i or --sequences_input SEQUENCES_INPUT
# File path to the raw sequences FASTA file to treat 
```

 ```shell
-db or --source_database {ncbi,ebi,ddbj}
# Name of the database from which the FASTA file was obtained.
 ```


 
  -fp FORWARD_PRIMER, --forward_primer FORWARD_PRIMER
                        File path to the forward primer FASTA file.
  -fmt FW_MISMATCH_TOL, --fw_mismatch_tol FW_MISMATCH_TOL
                        Number of mismatch which are accepted for the forward primer annealing.
  -rp REVERSE_PRIMER, --reverse_primer REVERSE_PRIMER
                        File path to the reverse primer FASTA file.
  -rmt RV_MISMATCH_TOL, --rv_mismatch_tol RV_MISMATCH_TOL
                        Number of mismatch which are accepted for the reverse primer annealing.



## Options

optional arguments:
  -h, --help            show this help message and exit
  -i SEQUENCES_INPUT, --sequences_input SEQUENCES_INPUT
                        File path to the raw sequences FASTA file to treat.
  -db {ncbi,ebi,ddbj}, --source_database {ncbi,ebi,ddbj}
                        Name of the database from which the FASTA file was obtained.
  -fp FORWARD_PRIMER, --forward_primer FORWARD_PRIMER
                        File path to the forward primer FASTA file.
  -fmt FW_MISMATCH_TOL, --fw_mismatch_tol FW_MISMATCH_TOL
                        Number of mismatch which are accepted for the forward primer annealing.
  -rp REVERSE_PRIMER, --reverse_primer REVERSE_PRIMER
                        File path to the reverse primer FASTA file.
  -rmt RV_MISMATCH_TOL, --rv_mismatch_tol RV_MISMATCH_TOL
                        Number of mismatch which are accepted for the reverse primer annealing.
  -inf INFOS_FILE, --infos_file INFOS_FILE
                        Informations file in which all the files treatments are recorded.
  -dit, --displ_inf_terminal
                        Display information file in terminal if arg is specified.
  -g1 GENUS1, --genus1 GENUS1
                        Name of the first genus to keep. Ex: Fusarium
  -g2 GENUS2, --genus2 GENUS2
                        Name of the second genus to keep. Ex: Gibberella
  -u, --unverified      Keep the unverified genus if arg is specified.
  -tbc, --taxonomy_by_complex
                        Do not group species with same complex lineage if arg is specified.
  -f {taxon,global}, --filtering_type {taxon,global}
  -rsv, --redund_seq_variants
                        Keep the redundant sequences variants if arg is specified.
  -ra, --redundant_amplicon
                        Keep the redundant amplicons if arg is specified.
  -nasv, --seq_without_ampl
                        Keep the sequences without amplicons if arg is specified.
  -svo SEQ_VARIANTS_OUTPUT, --seq_variants_output SEQ_VARIANTS_OUTPUT
                        File name of the sequences variants FASTA output.
  -ao AMPLICONS_OUTPUT, --amplicons_output AMPLICONS_OUTPUT
                        File name of the amplicons FASTA output.
  -to TAXONOMY_OUTPUT, --taxonomy_output TAXONOMY_OUTPUT
                        File name of the taxonomy TXT output.
  -tlo TAXON_LIST_OUTPUT, --taxon_list_output TAXON_LIST_OUTPUT
                        File name of the taxon list TXT output.
  -sao SHARED_AMPL_OUTPUT, --shared_ampl_output SHARED_AMPL_OUTPUT
                        File name of the shared amplicons TXT output.
  -cdo COMPLEX_DICT_OUTPUT, --complex_dict_output COMPLEX_DICT_OUTPUT
                        File name of the complex dictionnary TXT output.
  -mto MODIFIED_TAX_OUTPUT, --modified_tax_output MODIFIED_TAX_OUTPUT
                        File name of the modified taxon list TXT output.
  -ado ACCESS_DICT_OUTPUT, --access_dict_output ACCESS_DICT_OUTPUT
                        File name of the access dictionnary JSON output.

### 





# Importation pour module





