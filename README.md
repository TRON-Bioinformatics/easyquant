# EasyQuant (bp-quant)

*Quantification of reads at defined positions to verify custom input sequences.*

EasyQuant takes target sequence with defined positions or regions as input (e.g., breakpoints, splice junctions, retained introns) and quantifies supporting RNA-seq reads by reporting the number of overlapping 
reads (junction reads) and read pairs spanning the position or region (spanning pairs).

## Workflow

- Input:
    - Target sequences and positions of interest (CSV/TSV format)
    - FASTQ files or BAM file
- Convert target sequences to FASTA format (`bp_quant csv2fasta`)
- Map reads against sequences using STAR/Bowtie2
    - Generate index of sequences as reference (`bp_quant index`)
    - Map reads (`bp_quant align`)
- Count supporting reads using `bp_quant count`
- Output: 
    - Table with read counts per input sequence


 
## Installation

### Dependencies

Python and packages defined in the conda [environment.yml](https://github.com/TRON-Bioinformatics/easyquant/blob/master/environment.yml)

```
# Build conda environment or add the dependencies to your path
conda env create -f environment.yml -n easyquant_env
conda activate easyquant_env
```


### Install from PyPi

```
pip install bp-quant
```

### Install from GitHub

```
git clone https://github.com/TRON-Bioinformatics/easyquant.git
cd easyquant

python -m build
pip install dist/*.whl
```

## Usage

```
usage: bp_quant pipeline [-h] [-1 FQ1] [-2 FQ2] [-b BAM] -s SEQ_TAB -o OUTPUT_FOLDER [-d BP_DISTANCE] [--allow_mismatches] [--interval_mode] [--skip_singleton]
                         [-m {star,bowtie2}] [-t NUM_THREADS] [--alignment_params ALIGN_PARAMS] [--keep_aln | --keep_all]

Runs the complete bpquant pipeline

optional arguments:
  -h, --help            show this help message and exit
  -1 FQ1, --fq1 FQ1     Specify path to Read 1 (R1) FASTQ file
  -2 FQ2, --fq2 FQ2     Specify path to Read 2 (R2) FASTQ file
  -b BAM, --bam_file BAM
                        Specify path to input BAM file as alternative to FASTQ input
  -s SEQ_TAB, --sequence_tab SEQ_TAB
                        Specify the reference sequences as table with colums name, sequence, and position
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Specify the folder to save the results into.
  -d BP_DISTANCE, --bp_distance BP_DISTANCE
                        Threshold in base pairs for the required overlap size of reads on both sides of the breakpoint for junction/spanning read counting
  --allow_mismatches    Allow mismatches within the region around the breakpoint determined by the bp_distance parameter
  --interval_mode       Specify if interval mode shall be used
  --skip_singleton      Skip singleton alignments in requantification
  -m {star,bowtie2}, --method {star,bowtie2}
                        Specify alignment software to generate the index
  -t NUM_THREADS, --threads NUM_THREADS
                        Specify number of threads to use for the alignment
  --alignment_params ALIGN_PARAMS
                        Specify custom commandline parameters to use for the alignment
  --keep_aln            Do not delete alignment files during clean up step
  --keep_all            Do not perform clean up step after re-quantification

Copyright (c) 2024 TRON gGmbH (See LICENSE for licensing details)
```

**Note: For the quantification of splice junction sequences, we recommend performing the targeted alignment with strict parameters.**  

For bowtie2 we recommend the following additional alignment parameters: `--alignment_params "--dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.01"`

For STAR we recommend the following additional alignment parameters: `--alignment_params "--outFilterMismatchNoverReadLmax 0.3 --scoreDelOpen -2 --scoreInsOpen -2 --scoreDelBase -2 --scoreInsBase -2"`

### Use case with example data

Here, we use toy example data from the folder `example_data`. It consists of a table 
with input sequences and positions, as well as two fastq files / one BAM file. 

Example run with fastq files as input:

```
bp_quant pipeline \
  -1 example_data/example_rna-seq_R1_001.fastq.gz \
  -2 example_data/example_rna-seq_R1_001.fastq.gz \
  -s example_data/CLDN18_Context_seq.tsv \
  -d 10 \
  -o example_out \
  -m star \
  -t 6
```

Example run with BAM as input:

```
bp_quant pipeline \
  -b example_data/example_rna-seq.bam \
  -s example_data/CLDN18_Context_seq.tsv \
  -d 10 \
  -o example_out \
  -m star \
  -t 6
```

## Input

### Table with input sequences

The input target sequences should be provided as a Tab or `;`-separated table with
unique names and the relative (1-based) position of interest 
(breakpoint, junction) or interval of interest.

Example of an input table:

|name     | sequence      | position  |
|:--------|:--------------|:----------|
|seq1     | AACCGCCACCG   |5          |
|seq2     | GTCCGTTGGCG   |5          |
|seq3     | AACCGCCCTGT   |5          |
|seq4     | CGGCATCATCG   |0,5,10     |


#### Fastq files / BAM file

The sequencing data needs to be provided as paired-end fastq files or an unsorted BAM file (no multimappers). 

## Output format

The main output consists of the file `<OUTPUT_FOLDER>/quantification.tsv`. The table contains raw read counts for each input sequence.
The output folder contains additional files, such as  a table with mapping information of each mapped read (`read_info.tsv.gz`)

### Columns in output file `quantification.tsv`

#### Without interval mode

 - **name** name of the input sequence
 - **pos** position of interest relative to input sequence 
 - **junc** reads overlapping the position of interest
 - **span** read pairs spanning the position of interest
 - **anch** maximal number of bases next to position of interest that are overlaped by a single read
 - **a** reads mapping to sequence left of the position of interest
 - **b** reads mapping to sequence right of the position of interest

#### With interval mode

 - **name**   name of the input sequence
 - **interval** interval of interest relative to input sequence
 - **overlap_interval_end_reads** reads overlapping the end of the interval by at least `BP_DISTANCE` bases
 - **span_interval_end_pairs** read pairs spanning the end of the interval
 - **within_interval** reads mapping fully onto the interval
 - **coverage_perc** percentual coverage of the interval by aligned reads
 - **coverage_mean** average coverage per base for the interval (fold coverage)
 - **coverage_median** median coverage per base for the interval

### Example output

The output of the example run `<OUTPUT_FOLDER>/quantification.tsv` using a mismatch ratio of 0.05 (default) should look like this:


| name          | pos | junc | span | anch | a    | b    |
|---------------|-----|------|------|------|------|------|
| CLDN18_1      | 400 | 670  | 0    | 0    | 4140 | 10994|
| CLDN18_2      | 361 | 32   | 0    | 0    | 36   | 10994|
| CLDN18_total  | 400 | 612  | 2    | 0    | 11334| 14820|
| CLDN18_1_fake | 400 | 0    | 0    | 0    | 4132 | 14818|
| CLDN18_2_fake | 361 | 0    | 0    | 0    | 76   | 14818|
| HPRT1         | 400 | 76   | 0    | 0    | 1088 | 1332 |
| HPRT1_dup     | 400 | 76   | 0    | 0    | 1088 | 1332 |
| HPRT1_similar | 400 | 76   | 0    | 0    | 1064 | 814  |


Using the interval mode (`--interval_mode`) the output will look like the following:

| name          | interval | overlap_interval_end_reads | span_interval_end_pairs | within_interval | coverage_perc | coverage_mean | coverage_median |
|---------------|----------|----------------------------|-------------------------|-----------------|---------------|---------------|-----------------|
| CLDN18_1      | 0_400    | 670                        | 0                       | 4140            | 0.0           | 0.0           | 0.0             |
| CLDN18_1      | 400_786  | 0                          | 0                       | 10994           | 0.0           | 0.0           | 0.0             |
| CLDN18_2      | 0_361    | 32                         | 0                       | 36              | 0.0           | 0.0           | 0.0             |
| CLDN18_2      | 361_747  | 0                          | 0                       | 10994           | 0.0           | 0.0           | 0.0             |
| CLDN18_total  | 0_400    | 612                        | 2                       | 11334           | 0.0           | 0.0           | 0.0             |
| CLDN18_total  | 400_786  | 0                          | 0                       | 14820           | 0.0           | 0.0           | 0.0             |
| CLDN18_1_fake | 0_400    | 0                          | 0                       | 4132            | 0.0           | 0.0           | 0.0             |
| CLDN18_1_fake | 400_786  | 0                          | 0                       | 14818           | 0.0           | 0.0           | 0.0             |
| CLDN18_2_fake | 0_361    | 0                          | 0                       | 76              | 0.0           | 0.0           | 0.0             |
| CLDN18_2_fake | 361_747  | 0                          | 0                       | 14818           | 0.0           | 0.0           | 0.0             |
| HPRT1         | 0_400    | 76                         | 0                       | 1088            | 0.0           | 0.0           | 0.0             |
| HPRT1         | 400_793  | 0                          | 0                       | 1332            | 0.0           | 0.0           | 0.0             |
| HPRT1_dup     | 0_400    | 76                         | 0                       | 1088            | 0.0           | 0.0           | 0.0             |
| HPRT1_dup     | 400_793  | 0                          | 0                       | 1332            | 0.0           | 0.0           | 0.0             |
| HPRT1_similar | 0_400    | 76                         | 0                       | 1064            | 0.0           | 0.0           | 0.0             |
| HPRT1_similar | 400_793  | 0                          | 0                       | 814             | 0.0           | 0.0           | 0.0             |

Hint: This is just an example to illustrate the design of the table. The actual results may differ.


## Things to consider

EasyQuant supports two aligners, which have several differences:

- STAR: 
  - end-to-end alignment with no soft-clipping
  - Slow for many or large reference sequences
  - several parameters to optimize alignments, which can be customized with `--alignment_params`
- bowtie2: 
  - end-to-end alignment might lead to insertions where the context sequence starts/ends
  - faster than STAR for short reference sequences (index creation parameters are calculated automatically)


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## References

If you use EasyQuant in your research, please cite the following publication:

Franziska Lang, Patrick Sorn, Martin Suchan, Alina Henrich, Christian Albrecht, Nina Köhl, Aline Beicht, Pablo Riesgo-Ferreiro, Christoph Holtsträter, Barbara Schrörs, David Weber, Martin Löwer, Ugur Sahin, Jonas Ibn-Salem, **Prediction of tumor-specific splicing from somatic mutations as a source of neoantigen candidates**, Bioinformatics Advances, Volume 4, Issue 1, 2024, vbae080, https://doi.org/10.1093/bioadv/vbae080

```
@article{10.1093/bioadv/vbae080,
    author = {Lang, Franziska and Sorn, Patrick and Suchan, Martin and Henrich, Alina and Albrecht, Christian and Köhl, Nina and Beicht, Aline and Riesgo-Ferreiro, Pablo and Holtsträter, Christoph and Schrörs, Barbara and Weber, David and Löwer, Martin and Sahin, Ugur and Ibn-Salem, Jonas},
    title = "{Prediction of tumor-specific splicing from somatic mutations as a source of neoantigen candidates}",
    journal = {Bioinformatics Advances},
    volume = {4},
    number = {1},
    pages = {vbae080},
    year = {2024},
    month = {05},
    issn = {2635-0041},
    doi = {10.1093/bioadv/vbae080},
    url = {https://doi.org/10.1093/bioadv/vbae080},
    eprint = {https://academic.oup.com/bioinformaticsadvances/article-pdf/4/1/vbae080/58192195/vbae080.pdf},
}
```