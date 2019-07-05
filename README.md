# easyquant

*Quantification of reads at defined positions in custom input sequences.*


Given a gene fusion or splicing junction of interest, this tool can quantify RNA-seq reads 
supporting the breakpoint (junction) by quantifying reads that map to the breakpoint (junction reads) and
read pairs that span the breakpoint (spanning pairs).


## Workflow

- Input:
    - Tabel with target sequences and breakpoint position
    - fastq files
- Map reads against sequences using STARR
    - Generate Index of sequences as reference
    - Map reads
- Generate BED file with defined regions (left part and right part)
- Count reads using `requantify.py`
- Output: 
    - Table with read counts per input sequence and fastq sample

## Dependencies

 - Python 2.7.15
 - STARR
 - samtools (>= 1.9)
 
## Installation

```
git clone git@gitlab.rlp.net:tron/easyquant.git
```

## Usage example


```
python easy_quant.py \
  -i rna-seq_R1.fq.gz rna-seq_R1.fq.gz \
  -f sequence_table.tsv \
  -o output_folder
```


### Input

#### Table with input sequences

As input a tab-separated table should be given holding the target sequence 
with unique names and the relativ position of the brekapoint or fusion junction.

Example format of the input table:

|name     | sequence      | position  |
|:--------|:--------------|:----------|
|seq1     | AACCGCCACCG   |5          |
|seq2     | GTCCGTTGGCG   |5          |
|seq3     | AACCGCCCTGT   |5          |


#### Fastq files

Paired fastq files from the same sample should have an identical prefix in file 
name and contain `R1` or `R2` for the forward and reverse read pair, 
respectively. 

#### Config file


