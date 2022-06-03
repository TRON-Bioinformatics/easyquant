#!/bin/bash

#================================================================
# Test run for easy_qant.py on example data
#================================================================

# Remove existing output folder
rm -rf example_out_csv

# Run pipeline
python easy_quant.py \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out_csv

#================================================================
# Test run using a tab-separated input file
#================================================================

# Remove existing output folder
rm -rf example_out_tab

# Run pipeline
python easy_quant.py \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.tsv \
  -o example_out_tab


#================================================================
# Test run using different aligner
#================================================================

# Remove existing output folder
rm -rf example_out_csv_bowtie2

# Run pipeline
python easy_quant.py \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.tsv \
  -o example_out_csv_bowtie2 \
  -m bowtie2

# Remove existing output folder
rm -rf example_out_csv_bwa

# Run pipeline
python easy_quant.py \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.tsv \
  -o example_out_csv_bwa \
  -m bwa

#================================================================
# Test run using interval mode
#================================================================

# Remove existing output folder
rm -rf example_out_csv_interval

# Run pipeline
python easy_quant.py \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out_csv_interval \
  --interval-mode
