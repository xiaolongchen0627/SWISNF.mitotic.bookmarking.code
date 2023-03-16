bam=$1
PREFIX=$2
thread=$3
rsem-calculate-expression \
     --num-threads $thread \
     --no-bam-output  \
     --alignments  \
     --paired-end  \
     --strandedness reverse \
     $bam \
     RSEM-index/v1.3.1/GRCh38.primary_assembly.genome \
     ${PREFIX}.RSEM

getRSEMGeneCount.py \
     -g gencode.v37.primary_assembly.annotation.gtf.gene \
     -r ${PREFIX}