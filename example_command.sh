#!/bin/sh

python /home/sorn/development/ngs_pipelines/neo_orf_analysis/get_neo_orfs.py \
  -i /miseq/190313_M04039_0060_000000000-CB57P/Data/Intensities/BaseCalls/Library-1_S1_L001_R1_001.fastq.gz /miseq/190313_M04039_0060_000000000-CB57P/Data/Intensities/BaseCalls/Library-1_S1_L001_R2_001.fastq.gz \
  -f /scratch/info/projects/CM29_RNA_Seq/fusion_gene_detection/WP4_MiSeq/Sample_0207_003_BP_MET_1TRP_SUB13/Context_seqs.csv \
  -o /scratch/info/projects/CM29_RNA_Seq/fusion_gene_detection/WP4_MiSeq/Sample_0207_003_BP_MET_1TRP_SUB13/ \
  -p CentOS


# Innotop stomach samples
# IT_N_104
# IT_N_105
# and others

# Innotop Lung samples:
# IT_N_222
# IT_N_223
# and others

mkdir -p /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq

zcat /scratch/info/data/CM29_RNA_Seq/innotop_all/IT_N_222_SL1_S5_L005_R1_001.fastq.gz \
 | head -n 40000 \
 | gzip \
 > /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R1_001.h10k.fastq.gz

zcat /scratch/info/data/CM29_RNA_Seq/innotop_all/IT_N_222_SL1_S5_L005_R2_001.fastq.gz \
 | head -n 40000 \
 | gzip \
 > /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R2_001.h10k.fastq.gz


# build toy example fastqs fom mapping reads
module load bioinf/samtools/1.9
module load bioinf/bedtools/2.27.1
module load python/2.7.15

samtools sort -n /projects/CM29_RNA_Seq/easyquant/easyquant/Aligned.sortedByCoord.out.bam \
  -o /projects/CM29_RNA_Seq/easyquant/easyquant/Aligned.sortedByCoord.out.qsort.bam

bamToFastq -i /projects/CM29_RNA_Seq/easyquant/easyquant/Aligned.sortedByCoord.out.qsort.bam \
  -fq /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R1_001.mapping.fastq \
  -fq2 /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R2_001.mapping.fastq

gzip /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R1_001.mapping.fastq \
  > /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R1_001.mapping.fastq.gz
gzip /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R2_001.mapping.fastq \
  > /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R2_001.mapping.fastq.gz

# # samll referece genome
# python easy_quant.py \
#   -i /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R1_001.h10k.fastq.gz \
#      /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R2_001.h10k.fastq.gz \
#   -f /projects/CM29_RNA_Seq/easyquant/easyquant/example_data/CLDN18_Context_seqs_reformated.csv \
#   -o /projects/CM29_RNA_Seq/easyquant/easyquant/example_out_small

# samll data example
python easy_quant.py \
  -i /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R1_001.mapping.fastq.gz \
     /scratch/info/projects/CM29_RNA_Seq/easyquant/test_fq/IT_N_222_SL1_S5_L005_R2_001.mapping.fastq.gz \
  -f /projects/CM29_RNA_Seq/easyquant/easyquant/example_data/CLDN18_Context_seqs_reformated.csv \
  -o /projects/CM29_RNA_Seq/easyquant/easyquant/example_out_3

# 
# # full example
# python easy_quant.py \
#   -i /scratch/info/data/CM29_RNA_Seq/innotop_all/IT_N_222_SL1_S5_L005_R1_001.fastq.gz \
#      /scratch/info/data/CM29_RNA_Seq/innotop_all/IT_N_222_SL1_S5_L005_R2_001.fastq.gz \
#   -f /projects/CM29_RNA_Seq/easyquant/easyquant/example_data/CLDN18_Context_seqs_reformated_new.csv \
#   -o /projects/CM29_RNA_Seq/easyquant/easyquant/example_out_full
  
  
  