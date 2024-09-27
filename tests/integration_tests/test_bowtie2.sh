#!/bin/bash

#================================================================
# Test run using a tab-separated input file
# Input data: FASTQ
# Aligner: bowtie2
# Interval Mode: False
# Allow Mismatches in BP area(s): False
# Keep temporary files: False
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
