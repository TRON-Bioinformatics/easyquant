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
- Map reads against sequences using STAR
    - Generate Index of sequences as reference
    - Map reads
- Count reads using `requantify.py`
- Output: 
    - Table with read counts per input sequence and fastq sample

## Dependencies

 - Python 3
   - pysam (>= 0.16.0.1)
 - STAR (>= 2.6.1d)
 - samtools (>= 1.9)
 
## Installation

```
git clone https://github.com/TRON-Bioinformatics/easyquant.git
```

Update `config.ini` with installation paths

```
samtools_cmd=/path/to/samtools/1.9/samtools
star_cmd=/path/to/STAR/2.6.1d/bin/Linux_x86_64_static/STAR
```

## Usage


```
usage: easy_quant.py [-h] -i INPUT [INPUT ...] -s SEQ_TAB [-d BP_DISTANCE] -o
                     OUTPUT_FOLDER

Processing of demultiplexed FASTQs

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Specify the fastq folder(s) or fastq file(s) to
                        process.
  -s SEQ_TAB, --sequence_tab SEQ_TAB
                        Specify the reference sequences as table with colums
                        name, sequence, and position
  -d BP_DISTANCE, --bp_distance BP_DISTANCE
                        Threshold in base pairs for the required overlap size
                        of reads on both sides of the breakpoint for
                        junction/spanning read counting
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Specify the folder to save the results into.
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
  -o example_out
  
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

Paired fastq files from the same sample should have an identical prefix in file 
name and contain `R1` or `R2` for the forward and reverse read pair, 
respectively. 

#### Config file


### Output format

The main output consists of the file: 

 - `<OUTPUT_FOLDER>/quantification.tsv` contains raw read counts for each sequence

The output of the example data `<OUTPUT_FOLDER>/quantification.tsv` should look like this:


| name          | pos | junc | span | anch | a    | b    |
|:--------------|----:|-----:|-----:|-----:|-----:|-----:|
| CLDN18_1      | 400 | 570  | 689  | 25   | 1096 | 2497 |
| CLDN18_2      | 361 | 0    | 1    | 0    | 1    | 1529 |
| CLDN18_total  | 400 | 596  | 689  | 25   | 2770 | 3467 |
| CLDN18_1_fake | 400 | 2    | 3    | 14   | 5    | 2425 |
| CLDN18_2_fake | 361 | 0    | 0    | 0    | 0    | 2420 |
| HPRT1         | 400 | 107  | 216  | 25   | 848  | 686  |


Using the interval mode the output will look slightly different:

| name          | interval | overlap_stop | span_read | within_interval | coverage_perc | coverage_mean |
|:--------------|---------:|-------------:|----------:|----------------:|--------------:|--------------:|
| CLDN18_1      | 0_400    | 570          | 969       | 1191            | 0.89          | 191.9775      |
| CLDN18_1      | 400_786  | 0            | 969       | 3817            | 1.0           | 508.342       |
| CLDN18_2      | 0_361    | 0            | 1         | 1               | 0.141         | 0.141         |
| CLDN18_2      | 361_747  | 0            | 1         | 75              | 0.953         | 9.513         |
| CLDN18_total  | 0_400    | 596          | 1097      | 1624            | 1.0           | 244.4125      |
| CLDN18_total  | 400_786  | 0            | 1097      | 1409            | 1.0           | 183.847       |
| CLDN18_1_fake | 0_400    | 2            | 5         | 5               | 0.2525        | 0.705         |
| CLDN18_1_fake | 400_786  | 0            | 5         | 4199            | 1.0           | 548.430       |
| CLDN18_2_fake | 0_361    | 0            | 0         | 0               | 0.0           | 0.0           |
| CLDN18_2_fake | 361_747  | 0            | 0         | 356             | 1.0           | 41.003        |
| HPRT1         | 0_400    | 107          | 341       | 1443            | 1.0           | 187.12        |
| HPRT1         | 400_793  | 0            | 341       | 1082            | 1.0           | 138.483       |




#### Columns in output file

 - **name**   name of the input sequence
 - **pos** position of interest relative to input sequence 
 - **junc** reads overlapping the position of interest
 - **span** read pairs spanning the position of interest
 - **anch** maximal number of bases next to position of interest that are overlaped by a single read
 - **a** reads mapping to sequence left of the position of interest
 - **b** reads mapping to sequence right of the position of interest
 - **interval** interval of interest relative to input sequence
 - **overlap_stop** reads overlapping the end of the respective interval
 - **span_read** reads pairing with a read of a different interval
 - **within_interval** reads mapping fully onto the interval
 - **coverage_perc** percentual coverage of the interval by aligned reads
 - **coverage_mean** average coverage per base for the interval (fold coverage)