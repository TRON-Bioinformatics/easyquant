#!/bin/bash


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
  -d 10 \
  -o example_out_csv \
  -t 12 \
  -m star

#================================================================
# Test run using comma-separated input file
# Input data: FASTQ
# Aligner: STAR with individual parameters
# Interval Mode: False
# Allow Mismatches in BP area(s): False
# Keep temporary files: False
#================================================================

# Remove existing output folder
rm -rf example_out_csv_custom_params

# Run pipeline
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out_csv \
  -t 12 \
  -m star \
  --alignment_params "--outFilterMismatchNoverReadLmax 0.3 --scoreDelOpen -2 --scoreInsOpen -2 --scoreDelBase -2 --scoreInsBase -2"
