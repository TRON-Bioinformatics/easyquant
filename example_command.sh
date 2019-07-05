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


python easy_quant.py \
  -i /scratch/info/data/CM29_RNA_Seq/innotop_all/IT_N_222_SL1_S5_L005_R1_001.fastq.gz /scratch/info/data/CM29_RNA_Seq/innotop_all/IT_N_222_SL1_S5_L005_R2_001.fastq.gz \
  -f example_data/CLDN18_Context_seqs_reformated.csv \
  -o example_out