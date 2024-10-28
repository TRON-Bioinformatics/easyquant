#!/bin/bash

#================================================================
# Test run using a tab-separated input file
# Input data: FASTQ
# Aligner: STAR
# Interval Mode: False
# Allow Mismatches in BP area(s): False
# Keep temporary files: False
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
# Test run using comma-separated input file
# Input data: FASTQ
# Aligner: STAR
# Interval Mode: False
# Allow Mismatches in BP area(s): False
# Keep temporary files: False
#================================================================

# Remove existing output folder
rm -rf example_out_csv

# Run pipeline
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out_csv
