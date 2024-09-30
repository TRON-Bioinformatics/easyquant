#!/bin/bash

#================================================================
# Test run using a tab-separated input file
# Input data: FASTQ
# Aligner: bowtie2
# Interval Mode: False
# Allow Mismatches in BP area(s): False
# Keep temporary files: True
#================================================================


source tests/integration_tests/assert.sh
# Remove existing output folder
rm -rf example_out_csv_bowtie2
test_folder=example_out_csv_bowtie2
output=$test_folder/read_info.tsv.gz

# Run pipeline
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.tsv \
  -o $test_folder \
  -m bowtie2 \
  --keep_all

test -s $output || { echo "Missing read info file!"; exit 1; }
assert_eq `zcat $output | wc -l` 91919 "Wrong number of reads"

