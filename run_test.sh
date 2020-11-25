
module load  python/2.7.15

# Test run for requantify.py
# /code/Anaconda/3/2019/bin/python requantify.py \
python requantify.py \
  -i star_aln/Aligned.sortedByCoord.out.bam \
  -t example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o CLDN18_Context_seq.csv.counts_new5_star.tsv

# Test rund for easy_qant.py on example data
python easy_quant.py \
  -i example_data/example_rna-seq_R1_001.fastq.gz example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out_fix_insertion


# Run easyqant on example data
python easy_quant.py \
  -i example_data/example_rna-seq_R1_001.fastq.gz example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out_fix_insertion

python easy_quant.py \
  -i example_data/example_rna-seq_R1_001.fastq.gz example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out_0.2

# multiple sample input
python easy_quant.py \
  -i example_data/example_rna-seq_R1_001.fastq.gz \
      example_data/example_rna-seq_R2_001.fastq.gz \
      example_data/example_02_rna-seq_R1_001.fastq.gz \
      example_data/example_02_rna-seq_R2_001.fastq.gz \
  -d 3 \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out_mult

# mismatch rate example 
python easy_quant.py \
  -i example_data/example_rna-seq_R1_001.fastq.gz example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out_mmrate

