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
usage: easy_quant.py [-h] -1 FQ1 -2 FQ2 -s SEQ_TAB -o OUTPUT_FOLDER [-d BP_DISTANCE] [--allow_mismatches]
                     [--interval_mode] [-m {star,bowtie2,bwa}] [-t NUM_THREADS]

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
  --allow_mismatches    Allow mismatches within the region around the breakpoint determined by the bp_distance parameter
  --interval_mode       Specify if interval mode shall be used
  -m {star,bowtie2,bwa}, --method {star,bowtie2,bwa}
                        Specify alignment software to generate the index
  -t NUM_THREADS, --threads NUM_THREADS
                        Specify number of threads to use for the alignment


```

### Use case with example data

We use toy example data from the folder `example_data`. It consists of a table 
with input sequences and positions, as well as two fastq files. 

```
python easy_quant.py \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R2_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.csv \
  -d 10 \
  -o example_out \
  -m star \
  -t 6
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

The output of the example data `<OUTPUT_FOLDER>/quantification.tsv` using a mismatch ratio of 0.05 (default) should look like this:


| name          | pos | junc | span | anch | a    | b    |
|:--------------|----:|-----:|-----:|-----:|-----:|-----:|
| CLDN18_1      | 400 | 570  | 684  | 25   | 1109 | 3932 |
| CLDN18_2      | 361 | 0    | 1    | 0    | 1    | 2968 |
| CLDN18_total  | 400 | 593  | 688  | 25   | 4315 | 4990 |
| CLDN18_1_fake | 400 | 0    | 3    | 0    | 3    | 3946 |
| CLDN18_2_fake | 361 | 0    | 0    | 0    | 0    | 3943 |
| HPRT1         | 400 | 107  | 215  | 25   | 1259 | 974  |


Using the interval mode the output will look slightly different:

| name          | interval | overlap_interval_end_reads | span_interval_end_pairs | within_interval | coverage_perc | coverage_mean | coverage_median |
|:--------------|---------:|---------------------------:|------------------------:|----------------:|--------------:|--------------:|----------------:|
| CLDN18_1      | 0_400    | 570                        | 684                     | 1109            | 0.89          | 182.71        | 137.5           |
| CLDN18_1      | 400_786  | 0                          | 0                       | 3932            | 1.0           | 519.51        | 563.5           |
| CLDN18_2      | 0_361    | 0                          | 1                       | 1               | 0.14          | 0.14          | 0.0             |
| CLDN18_2      | 361_747  | 0                          | 0                       | 2968            | 1.0           | 392.15        | 425.5           |
| CLDN18_total  | 0_400    | 593                        | 688                     | 4315            | 1.0           | 586.16        | 682.0           |
| CLDN18_total  | 400_786  | 0                          | 0                       | 4990            | 1.0           | 659.15        | 757.5           |
| CLDN18_1_fake | 0_400    | 0                          | 3                       | 3               | 0.16          | 0.38          | 0.0             |
| CLDN18_1_fake | 400_786  | 0                          | 0                       | 3946            | 1.0           | 521.23        | 551.5           |
| CLDN18_2_fake | 0_361    | 0                          | 0                       | 0               | 0.0           | 0.0           | 0.0             |
| CLDN18_2_fake | 361_747  | 0                          | 0                       | 3943            | 1.0           | 520.83        | 551.0           |
| HPRT1         | 0_400    | 107                        | 215                     | 1259            | 1.0           | 167.77        | 175.0           |
| HPRT1         | 400_793  | 0                          | 0                       | 974             | 0.98          | 126.40        | 101.0           |




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