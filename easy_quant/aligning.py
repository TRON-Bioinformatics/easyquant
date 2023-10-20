#!/usr/bin/env python


import subprocess


def get_align_cmd_bowtie2(fq1, fq2, index_dir, out_sam, num_threads):
    cmd = "bowtie2 -p {0} -x {1}/bowtie -1 {2} -2 {3} -S {4}".format(
        num_threads,
        index_dir,
        fq1,
        fq2,
        out_sam
    )

    return cmd

    
def get_align_cmd_bwa(fq1, fq2, index_dir, out_sam, num_threads):
    cmd = "bwa mem -t {0} {1} {2} {3} > {4}".format(
        num_threads,
        index_dir, 
        fq1, 
        fq2, 
        out_sam
    )

    return cmd


def get_align_cmd_star(fq1, fq2, index_dir, out_sam, num_threads, custom_params):
    cmd = "STAR --outFileNamePrefix {0} \
    --limitOutSAMoneReadBytes 1000000 \
    --genomeDir {1} \
    --readFilesCommand 'gzip -d -c -f' \
    --readFilesIn {2} {3} \
    --outSAMmode Full \
    --alignEndsType EndToEnd \
    --outFilterMultimapNmax -1 \
    --outSAMattributes NH HI AS nM NM MD \
    --outSAMunmapped Within KeepPairs \
    --outFilterScoreMinOverLread 0.3 \
    --outFilterMatchNminOverLread 0.3 \
    {4} \
    --runThreadN {5}".format(
        out_sam, 
        index_dir, 
        fq1,
        fq2,
        custom_params,
        num_threads
    )


    return cmd


def run(fq1, fq2, index_dir, out_path, threads, method, params):
    if method == "bowtie2":
        subprocess.run(get_align_cmd_bowtie2(fq1, fq2, index_dir, out_path, threads).split(" "))
    elif method == "bwa":
        subprocess.run(get_align_cmd_bwa(fq1, fq2, index_dir, out_path, threads).split(" "))
    elif method == "star":
        subprocess.run(get_align_cmd_star(fq1, fq2, index_dir, out_path, threads, params).split(" "))


def add_aligner_args(parser):
    parser.add_argument(
        "-1",
        "--fq1",
        dest="fq1",
        help="Specify FQ1",
        required=True
    )
    parser.add_argument(
        "-2",
        "--fq2",
        dest="fq2",
        help="Specify FQ2",
        required=True
    )
    parser.add_argument(
        "-i",
        "--index_dir",
        dest="index_dir",
        help="Specify path to input index dir",
        required=True
    )
    parser.add_argument(
        "-o",
        "--output_path",
        dest="output_path",
        help="Specify path to output file/folder",
        required=True
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        help="Specify amount of threads to be used during runtime",
        default=1,
    )
    parser.add_argument(
        "-m",
        "--method",
        dest="method",
        choices=["star", "bowtie2", "bwa"],
        help="Specify aligner for execution",
        default="star",
    )
    parser.add_argument(
        "-p",
        "--params",
        dest="params",
        help="Specify custom params (e.g. for STAR)"
    )
    parser.set_defaults(func=alignment_command)


def alignment_command(args):
    run(args.fq1, args.fq2, args.index_dir, args.output_path, args.threads, args.method, args.params)
