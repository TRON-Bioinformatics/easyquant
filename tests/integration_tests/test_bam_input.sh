#!/bin/bash

#================================================================
# Test run using comma-separated input file
# Input data: BAM
# Aligner: STAR
# Interval Mode: False
# Allow Mismatches in BP area(s): True
# Keep temporary files: False
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
