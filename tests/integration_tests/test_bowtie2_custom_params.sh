#!/bin/bash


#================================================================
# Test run using a tab-separated input file
# Input data: FASTQ
# Aligner: bowtie2 with individual parameters
# Interval Mode: False
# Allow Mismatches in BP area(s): False
# Keep temporary files: False
#================================================================

# Remove existing output folder
rm -rf example_out_csv_bowtie2_1

# Run pipeline
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.tsv \
  -o example_out_csv_bowtie2_1 \
  -m bowtie2 \
  --alignment_params "--dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.01"



#================================================================
# Test run using a tab-separated input file
# Input data: FASTQ
# Aligner: bowtie2 with individual parameters
# Interval Mode: False
# Allow Mismatches in BP area(s): False
# Keep temporary files: False
#================================================================

# Remove existing output folder
rm -rf example_out_csv_bowtie2_2

# Run pipeline
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.tsv \
  -o example_out_csv_bowtie2_2 \
  -m bowtie2 \
  --alignment_params "--dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.02"
