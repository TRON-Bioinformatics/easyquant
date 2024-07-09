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


#================================================================
# Test run using a tab-separated input file
# Input data: FASTQ
# Aligner: bwa
# Interval Mode: False
# Allow Mismatches in BP area(s): False
# Keep temporary files: False
#================================================================

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
# Test run using comma-separated input file
# Input data: FASTQ
# Aligner: STAR
# Interval Mode: True
# Allow Mismatches in BP area(s): False
# Keep temporary files: False
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
# Test run using comma-separated input file
# Input data: FASTQ
# Aligner: STAR
# Interval Mode: False
# Allow Mismatches in BP area(s): True
# Keep temporary files: False
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
