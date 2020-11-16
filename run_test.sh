
module load  python/2.7.15

# Run example data 
python easy_quant.py \
  -i example_data/example_rna-seq_R1_001.fastq.gz example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out

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

