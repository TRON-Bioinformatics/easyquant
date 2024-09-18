import subprocess

def get_align_cmd_bowtie2(fq1, fq2, bam, index_dir, out_dir, num_threads, custom_params):
    cmd=""
    if fq1 and fq2:
        cmd = "bowtie2 -p {0} -x {1}/bowtie -a --end-to-end --no-discordant -1 {2} -2 {3} -S {4}/Aligned.out.sam {5}".format(
            num_threads,
            index_dir,
            fq1,
            fq2,
            out_dir,
            custom_params
        )
    elif not fq1 and not fq2 and bam:
        cmd = "bowtie2 -p {0} -x {1}/bowtie -a --end-to-end --no-discordant -b {2} --align-paired-reads -S {3}/Aligned.out.sam {4}".format(
            num_threads,
            index_dir,
            bam,
            out_dir,
            custom_params
        ) 
    return cmd

def get_align_cmd_star(fq1, fq2, bam, index_dir, out_dir, num_threads, custom_params):
    cmd = ""
    if fq1 and fq2:
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
            out_dir + "/",
            index_dir,
            fq1,
            fq2,
            custom_params,
            num_threads
        )
    elif not fq1 and not fq2 and bam:
        cmd = "STAR --outFileNamePrefix {0} \
        --runMode alignReads \
        --limitOutSAMoneReadBytes 1000000 \
        --genomeDir {1} \
        --readFilesType SAM PE \
        --readFilesCommand 'samtools view' \
        --readFilesIn {2} \
        --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
        --outSAMmode Full \
        --alignEndsType EndToEnd \
        --outFilterMultimapNmax -1 \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMunmapped Within KeepPairs \
        --outFilterScoreMinOverLread 0.3 \
        --outFilterMatchNminOverLread 0.3 \
        {3} \
        --runThreadN {4}".format(
            out_dir + "/",
            index_dir,
            bam,
            custom_params,
            num_threads
        )
    return cmd


def run(fq1, fq2, bam, index_dir, out_path, threads, method, params):
    cmd = None
    if method == "bowtie2":
        cmd = get_align_cmd_bowtie2(fq1, fq2, bam, index_dir, out_path, threads, params)
    elif method == "star":
        cmd = get_align_cmd_star(fq1, fq2, bam, index_dir, out_path, threads, params)
    print(cmd)
    subprocess.run(cmd, shell=True, cwd=out_path)

def add_aligner_args(parser):
    parser.add_argument(
        "-1",
        "--fq1",
        dest="fq1",
        help="Specify FQ1"
    )
    parser.add_argument(
        "-2",
        "--fq2",
        dest="fq2",
        help="Specify FQ2"
    )
    parser.add_argument(
        "-b",
        "--bam",
        dest="bam",
        help="Specify unaligned input BAM file"
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
        choices=["star", "bowtie2"],
        help="Specify aligner for execution",
        default="star",
    )
    parser.add_argument(
        "-p",
        "--params",
        dest="params",
        help="Specify custom params (e.g. for STAR)",
        default=""
    )
    parser.set_defaults(func=alignment_command)


def alignment_command(args):
    run(
        fq1=args.fq1,
        fq2=args.fq2,
        bam=args.bam,
        index_dir=args.index_dir,
        out_path=args.output_path,
        threads=args.threads,
        method=args.method,
        params=args.params
    )
