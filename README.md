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

 - Python
 - STARR (for read mapping)

## Installation

```
git clone git@gitlab.rlp.net:tron/easyquant.git
```

## Usage example


```
python <script> <args>
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

#### Config file


