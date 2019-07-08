# easyquant

*Quantification of reads at defined positions to varify custom input sequences.*


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

## Usage


```
python easy_quant.py -h
usage: easy_quant.py [-h] -i INPUT [INPUT ...] -s SEQ_TAB -o OUTPUT_FOLDER

Processing of demultiplexed FASTQs

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Specify the fastq folder(s) or fastq file(s) to
                        process.
  -s SEQ_TAB, --sequence_tab SEQ_TAB
                        Specify the reference seqeunces as table with colums
                        name, sequence, and position
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Specify the folder to save the results into.
```

### Use case with example data

We use toy example data from the folder `example_data`. It consists of a table 
with input sequences and positions, as well as two fastq files. 

```
python easy_quant.py \
  -i example_data/example_rna-seq_R1_001.fastq.gz example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -o example_out
  
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


### Output format

The main output consists of two files: 

 - `<OUTPUT_FOLDER>/quantification.csv.counts` contains normalized read counts for each sequence
 - `<OUTPUT_FOLDER>/quantification.csv.counts` contains raw read counts for each sequence

The output of the example data `<OUTPUT_FOLDER>/quantification.csv.counts` should look like this:

|name          | position|    a|    b| junc| span| anch|
|:-------------|--------:|----:|----:|----:|----:|----:|
|CLDN18_1      |      400| 1097| 2500|  685|  693|   25|
|HPRT1         |      400|  848|  690|  200|  220|   25|
|CLDN18_2      |      361|    1| 1529|   21|    1|    7|
|CLDN18_2_fake |      361|    0| 2420|   83|    0|    7|
|CLDN18_total  |      400| 2794| 3483|  834|  729|   25|
|CLDN18_1_fake |      400|    5| 2425|   87|    3|   14|


#### Columns in ouptut file

 - **name**   name of the input sequence
 - **position** position of interest relative to input sequence 
 - **a** reads mapping to sequence left of the position of interest
 - **b** reads mapping to sequence right of the position of interest
 - **junc** reads overalpping the position of interest
 - **span** read pairs spanning the position of interest
 - **anch** maximal number of bases next to position of interest that are overlaped by a single read

