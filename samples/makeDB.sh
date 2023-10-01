GENUS=$1
DB=$2
SEQ=$3

mkdir database

#=================================
#  1         CREATE
#=================================

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

#=================================
#  2         FILTER
#=================================
python3 -m asvmaker \
    -inf database/${DB}/${GENUS}_${DB}_info_filter.txt \
    filter \
    -i database/${DB}/${GENUS}_${DB}_create.json \
    -g1 ${GENUS} \
    -o database/${DB}/${GENUS}_${DB}_filter.json

#=================================
#  3         EXPORT 1
#=================================
python3 -m asvmaker \
    -inf database/${DB}/${GENUS}_${DB}_info_exp1.txt \
    export \
    -i database/${DB}/${GENUS}_${DB}_filter.json \
    -sao database/${DB}/${GENUS}_${DB}_sa.txt 1

#=================================
#  4         EDIT
#=================================
python3 -m asvmaker \
    -inf database/${DB}/${GENUS}_${DB}_info_edit.txt \
    edit \
    -i database/${DB}/${GENUS}_${DB}_filter.json \
    -grp database/${DB}/${GENUS}_${DB}_sa_ext.txt \
    -o database/${DB}/${GENUS}_${DB}_edit.json

#=================================
#  5         EXPORT 2
#=================================
python3 -m asvmaker \
    -inf database/${DB}/${GENUS}_${DB}_info_exp2.txt \
    export \
    -i database/${DB}/${GENUS}_${DB}_edit.json \
    -aop database/${DB}/${GENUS}_${DB}_asv.fasta