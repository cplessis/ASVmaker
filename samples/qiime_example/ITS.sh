#==============================================
#          Parameters
#==============================================

INPUT_FASTQ=$1
OUTPUT_PATH=$2
REF_PATH=$3

CLASSIFIER=UNITE8.3_classifier.qza
REF_DB=PathDB_BITS_2021.6-phylo.fasta
REF_DB_TAXO=Genus_Fungi_taxo.txt

echo "Input files : "  $INPUT_FASTQ
echo "Output path : "  $OUTPUT_PATH
echo "References path : "  $REF_PATH
echo "References database : "  $REF_DB
echo "References database taxo : "  $REF_DB_TAXO
echo "Classifier : "  $CLASSIFIER


#==============================================
#          Qiime Analysis
#==============================================


# Import FASTQ files
qiime tools import \
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path $INPUT_FASTQ \
 --input-format CasavaOneEightSingleLanePerSampleDirFmt \
 --output-path $OUTPUT_PATH/demux-paired-end.qza


# Trim primers
qiime cutadapt trim-paired \
    --p-cores 6 \
    --i-demultiplexed-sequences $OUTPUT_PATH/demux-paired-end.qza \
    --p-adapter-f '^ACCTGCGGARGGATCA...AACTTTYARCAAYGGATCTC' \
    --p-adapter-r '^GAGATCCRTTGYTRAAAGTT...TGATCCYTCCGCAGGT' \
    --p-error-rate 0.2 \
    --o-trimmed-sequences $OUTPUT_PATH/cuta-trimmed.qza \
    --p-discard-untrimmed \
    --p-match-read-wildcards

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $OUTPUT_PATH/cuta-trimmed.qza \
  --o-table $OUTPUT_PATH/table.qza \
  --o-representative-sequences $OUTPUT_PATH/rep-seqs-all.qza \
  --o-denoising-stats $OUTPUT_PATH/denoising-stats.qza \
  --p-n-threads 12 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0

qiime tools export \
  --input-path $OUTPUT_PATH/rep-seqs-all.qza \
  --output-path $OUTPUT_PATH

# Export to TSV
qiime tools export \
  --input-path $OUTPUT_PATH/table.qza \
  --output-path $OUTPUT_PATH

biom convert \
  -i $OUTPUT_PATH/feature-table.biom \
  -o $OUTPUT_PATH/table.tsv \
  --to-tsv

# Classifcation
qiime feature-classifier classify-sklearn \
  --i-classifier $REF_PATH/$CLASSIFIER \
  --i-reads $OUTPUT_PATH/rep-seqs-all.qza \
  --o-classification $OUTPUT_PATH/taxo-all.qza

#Exportation to TSV
qiime tools export \
  --input-path $OUTPUT_PATH/taxo-all.qza \
  --output-path $OUTPUT_PATH

#==============================================
#   ASV identification by ASVmaker specific database
#==============================================

python3 ./pathoPipeline/tools/qiime2csv.py \
    --outdir $OUTPUT_PATH \
    --kingdom k__Fungi \
    --dna-sequences $OUTPUT_PATH/dna-sequences.fasta \
    --asv-db-phylo $REF_PATH/$REF_DB \
    --asv-db-taxo $REF_PATH/$REF_DB_TAXO \
    --taxonomy $OUTPUT_PATH/taxonomy.tsv \
    --table $OUTPUT_PATH/table.tsv
