#!/bin/sh

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


# full example
python easy_quant.py \
  -i /scratch/info/data/CM29_RNA_Seq/innotop_all/IT_N_222_SL1_S5_L005_R1_001.fastq.gz \
     /scratch/info/data/CM29_RNA_Seq/innotop_all/IT_N_222_SL1_S5_L005_R2_001.fastq.gz \
  -s /projects/CM29_RNA_Seq/easyquant/easyquant/example_data/CLDN18_Context_seqs_reformated_new.csv \
  -o /projects/CM29_RNA_Seq/easyquant/easyquant/example_out_full
  
  

samtools sort -n example_out_full/star/Aligned.sortedByCoord.out.bam \
  -o example_out_full/star/Aligned.sortedByCoord.out.qsort.bam

bamToFastq -i example_out_full/star/Aligned.sortedByCoord.out.qsort.bam \
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
  -s /projects/CM29_RNA_Seq/easyquant/easyquant/example_data/example_data/CLDN18_Context_seq.csv \
  -o /projects/CM29_RNA_Seq/easyquant/easyquant/example_out_4


  