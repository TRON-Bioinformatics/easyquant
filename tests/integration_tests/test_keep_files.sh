#!/bin/bash

#================================================================
# Test run using a comma-separated input file
# Input data: FASTQ
# Aligner: STAR with individual parameters
# Interval Mode: False
# Allow Mismatches in BP area(s): False
# Keep temporary files: All
#================================================================

# Remove existing output folder
rm -rf example_out_csv_keep_all

# Run pipeline
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out_csv_keep_all \
  -t 12 \
  --keep_all


#================================================================
# Test run using comma-separated input file
# Input data: BAM
# Aligner: STAR
# Interval Mode: False
# Allow Mismatches in BP area(s): True
# Keep temporary files: CRAM only
#================================================================

# Remove existing output folder
rm -rf example_out_csv_keep_aln

# Run pipeline
bp_quant pipeline \
  -b example_data/example_rna-seq.bam \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out_csv_keep_aln \
  -t 12 \
  --keep_aln
