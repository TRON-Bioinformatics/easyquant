
# module load anaconda/3/2019

#================================================================
# Test run for easy_qant.py on example data
#================================================================
python easy_quant.py \
  -i example_data/example_rna-seq_R1_001.fastq.gz example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out_csv

#================================================================
# Test run using a tab-separated input file
#================================================================
python easy_quant.py \
  -i example_data/example_rna-seq_R1_001.fastq.gz example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.tsv \
  -o example_out_tab
