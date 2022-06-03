# easyquant

*Quantification of reads at defined positions to verify custom input sequences.*

Given a gene fusion or splicing junction of interest, this tool can quantify
RNA-seq reads supporting the breakpoint (or splice junction) by quantifying
reads that map to the breakpoint (junction reads) and read pairs that span the
breakpoint (spanning pairs).


## Workflow

- Input:
    - Table with target sequences and breakpoints position (CSV/TSV format)
    - fastq files
- Map reads against sequences using STAR/Bowtie2/BWA
    - Generate Index of sequences as reference
    - Map reads
- Count reads using `requantify.py`
- Output: 
    - Table with read counts per input sequence

## Dependencies

 - Python 3
   - pysam (>= 0.16.0.1)
 - STAR (>= 2.6.1d)
 - samtools (>= 1.9)
 
## Installation

```
git clone https://github.com/TRON-Bioinformatics/easyquant.git
```

```
mv config.ini.sample config.ini
```

Update `config.ini` with installation paths

```
samtools_cmd=/path/to/samtools/1.9/samtools
star_cmd=/path/to/STAR/2.6.1d/bin/Linux_x86_64_static/STAR
```

## Usage


```
usage: easy_quant.py [-h] -1 FQ1 -2 FQ2 -s SEQ_TAB -o OUTPUT_FOLDER [-d BP_DISTANCE] [-m {star,bowtie2,bwa}]
                     [--interval-mode]

Processing of demultiplexed FASTQs

optional arguments:
  -h, --help            show this help message and exit
  -1 FQ1, --fq1 FQ1     Specify path to Read 1 (R1) FASTQ file
  -2 FQ2, --fq2 FQ2     Specify path to Read 2 (R2) FASTQ file
  -s SEQ_TAB, --sequence_tab SEQ_TAB
                        Specify the reference sequences as table with colums name, sequence, and position
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Specify the folder to save the results into.
  -d BP_DISTANCE, --bp_distance BP_DISTANCE
                        Threshold in base pairs for the required overlap size of reads on both sides of the breakpoint for
                        junction/spanning read counting
  -m {star,bowtie2,bwa}, --method {star,bowtie2,bwa}
                        Specify alignment software to generate the index
  --interval-mode       Specify if interval mode shall be used

```

### Use case with example data

We use toy example data from the folder `example_data`. It consists of a table 
with input sequences and positions, as well as two fastq files. 

```
python easy_quant.py \
  -i example_data/example_rna-seq_R1_001.fastq.gz example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out \
  -m star
  [--interval-mode]
```



### Input

#### Table with input sequences

As input a CSV/TSV table should be given holding the target sequence 
with unique names and the relative position(s) of the breakpoint(s)/intervals or fusion junction.

Example format of the input table:

|name     | sequence      | position  |
|:--------|:--------------|:----------|
|seq1     | AACCGCCACCG   |5          |
|seq2     | GTCCGTTGGCG   |5          |
|seq3     | AACCGCCCTGT   |5          |
|seq4     | CGGCATCATCG   |0,5,10     |


#### Fastq files

Paired fastq files as input is required 
to successfully determine read classes as described below. 

#### Config file


### Output format

The main output consists of the file: 

 - `<OUTPUT_FOLDER>/quantification.tsv` contains raw read counts for each sequence

The output of the example data `<OUTPUT_FOLDER>/quantification.tsv` should look like this:


| name          | pos | junc | span | anch | a    | b    |
|:--------------|----:|-----:|-----:|-----:|-----:|-----:|
| CLDN18_1      | 400 | 570  | 689  | 25   | 1116 | 3989 |
| CLDN18_2      | 361 | 0    | 1    | 0    | 1    | 3021 |
| CLDN18_total  | 400 | 596  | 689  | 25   | 4373 | 5803 |
| CLDN18_1_fake | 400 | 2    | 3    | 14   | 5    | 4761 |
| CLDN18_2_fake | 361 | 0    | 0    | 0    | 0    | 4756 |
| HPRT1         | 400 | 107  | 216  | 25   | 1400 | 1021 |


Using the interval mode the output will look slightly different:

| name          | interval | overlap_interval_end_reads | span_interval_end_pairs | within_interval | coverage_perc | coverage_mean | coverage_median |
|:--------------|---------:|---------------------------:|------------------------:|----------------:|--------------:|--------------:|----------------:|
| CLDN18_1      | 0_400    | 570                        | 689                     | 1116            | 0.89          | 183.38        | 137.5           |
| CLDN18_1      | 400_786  | 0                          | 0                       | 3989            | 1.0           | 530.74        | 563.5           |
| CLDN18_2      | 0_361    | 0                          | 1                       | 1               | 0.14          | 0.14          | 0.0             |
| CLDN18_2      | 361_747  | 0                          | 0                       | 3021            | 1.0           | 402.0         | 425.5           |
| CLDN18_total  | 0_400    | 596                        | 689                     | 4373            | 1.0           | 599.54        | 682.0           |
| CLDN18_total  | 400_786  | 0                          | 0                       | 5803            | 1.0           | 754.26        | 757.5           |
| CLDN18_1_fake | 0_400    | 2                          | 3                       | 5               | 0.25          | 0.71          | 0.0             |
| CLDN18_1_fake | 400_786  | 0                          | 0                       | 4761            | 1.0           | 616.93        | 551.5           |
| CLDN18_2_fake | 0_361    | 0                          | 0                       | 0               | 0.0           | 0.0           | 0.0             |
| CLDN18_2_fake | 361_747  | 0                          | 0                       | 4756            | 1.0           | 616.27        | 551.0           |
| HPRT1         | 0_400    | 107                        | 216                     | 1400            | 1.0           | 182.21        | 175.0           |
| HPRT1         | 400_793  | 0                          | 0                       | 1021            | 1.0           | 131.49        | 101.0           |




#### Columns in output file

### While not using interval mode

 - **name**   name of the input sequence
 - **pos** position of interest relative to input sequence 
 - **junc** reads overlapping the position of interest
 - **span** read pairs spanning the position of interest
 - **anch** maximal number of bases next to position of interest that are overlaped by a single read
 - **a** reads mapping to sequence left of the position of interest
 - **b** reads mapping to sequence right of the position of interest

### While using interval mode

 - **name**   name of the input sequence
 - **interval** interval of interest relative to input sequence
 - **overlap_interval_end_reads** reads overlapping the end of the interval by at least `BP_DISTANCE` bases
 - **span_interval_end_pairs** read pairs spanning the end of the interval
 - **within_interval** reads mapping fully onto the interval
 - **coverage_perc** percentual coverage of the interval by aligned reads
 - **coverage_mean** average coverage per base for the interval (fold coverage)
 - **coverage_median** median coverage per base for the interval