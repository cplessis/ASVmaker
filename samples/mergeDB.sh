GENUS=$1

mkdir database/all

#=================================
#  1         MERGE
#=================================
python3 -m asvmaker \
    -inf database/all/${GENUS}_info_merge.txt \
    merge \
    -i database/rnaCentral/${GENUS}_rnaCentral_edit.json \
    -sa1 database/rnaCentral/${GENUS}_rnaCentral_sa_ext.txt \
    -i2 database/unite/${GENUS}_unite_edit.json \
    -sa2 database/unite/${GENUS}_unite_sa_ext.txt \
    -o database/all/${GENUS}_merge.json

#=================================
#  2         EXPORT
#=================================
python3 -m asvmaker \
    -inf database/all/${GENUS}_info_exp3.txt \
    export \
    -i database/all/${GENUS}_merge.json \
    -aop database/all/${GENUS}_asv.fasta

    