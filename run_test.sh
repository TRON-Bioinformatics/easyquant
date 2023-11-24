#!/bin/bash

#================================================================
# Test run for bp_qant.py on example data
#================================================================

# Remove existing output folder
rm -rf example_out_csv

# Run pipeline
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out_csv \
  -t 12 \
  --star_cmd_params "--outFilterMismatchNoverReadLmax 0.3 --scoreDelOpen -2 --scoreInsOpen -2 --scoreDelBase -2 --scoreInsBase -2"

#================================================================
# Test run using a tab-separated input file
#================================================================

# Remove existing output folder
rm -rf example_out_tab

# Run pipeline
bp_quant pipeline \
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
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.tsv \
  -o example_out_csv_bowtie2 \
  -m bowtie2

# Remove existing output folder
rm -rf example_out_csv_bwa

# Run pipeline
bp_quant pipeline \
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
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out_csv_interval \
  --interval_mode


#================================================================
# Test run using interval mode
#================================================================

# Remove existing output folder
rm -rf example_out_csv_mismatch

# Run pipeline
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out_csv_mismatch \
  --allow_mismatches


#================================================================
# Test run using BAM file as input
#================================================================

# Remove existing output folder
rm -rf example_out_csv_bam

# Run pipeline
bp_quant pipeline \
  -b example_data/example_rna-seq.bam \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out_csv_bam \
  -t 12


#================================================================
# Test run keeping alignment as BAM
#================================================================

# Remove existing output folder
rm -rf example_out_keep_bam

# Run pipeline
bp_quant pipeline \
  -b example_data/example_rna-seq.bam \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out_keep_bam \
  -t 12 \
  --keep_bam

#================================================================
# Test run without clean up step
#================================================================

# Remove existing output folder
rm -rf example_out_keep_all

# Run pipeline
bp_quant pipeline \
  -b example_data/example_rna-seq.bam \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out_keep_all \
  -t 12 \
  --keep_all
