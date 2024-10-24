"""
This module generates the alignment commands depending on the input parameters.
"""

# pylint: disable=E0401
import bp_quant.io.io_methods as IOMethods

def get_align_cmd_bowtie2(fq1, fq2, bam, index_dir, out_dir, num_threads, custom_params) -> str:
    """Build up alignment command for bowtie2.

    Args:
        fq1 (str): Input FQ1 file
        fq2 (str): Input FQ2 file
        bam (str): Input BAM file
        index_dir (str): Index path to use for processing
        out_dir (str): Output directory where results are stored to
        num_threads (int): Number of cpus to use for processing
        custom_params (str): Custom parameters string to add to the final command

    Returns:
        str: command to be executed in the subprocess
    """
    cmd=""
    if fq1 and fq2:
        cmd = f"bowtie2 \
            -p {num_threads} \
            -x {index_dir}/bowtie \
            -a --end-to-end --no-discordant \
            -1 {fq1} -2 {fq2} \
            -S {out_dir}/Aligned.out.sam \
            {custom_params}"
    elif not fq1 and not fq2 and bam:
        cmd = f"bowtie2 \
            -p {num_threads} \
            -x {index_dir}/bowtie \
            -a --end-to-end --no-discordant \
            -b {bam} \
            --align-paired-reads \
            -S {out_dir}/Aligned.out.sam \
            {custom_params}"
    return cmd

def get_align_cmd_star(fq1, fq2, bam, index_dir, out_dir, num_threads, custom_params) -> str:
    """Build up alignment command for STAR.

    Args:
        fq1 (str): Input FQ1 file
        fq2 (str): Input FQ2 file
        bam (str): Input BAM file
        index_dir (str): Index path to use for processing
        out_dir (str): Output directory where results are stored to
        num_threads (int): Number of cpus to use for processing
        custom_params (str): Custom parameters string to add to the final command

    Returns:
        str: command to be executed in the subprocess
    """
    cmd = ""
    if fq1 and fq2:
        cmd = f"STAR --outFileNamePrefix {out_dir}/ \
        --limitOutSAMoneReadBytes 1000000 \
        --genomeDir {index_dir} \
        --readFilesCommand 'gzip -d -c -f' \
        --readFilesIn {fq1} {fq2} \
        --outSAMmode Full \
        --alignEndsType EndToEnd \
        --outFilterMultimapNmax -1 \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMunmapped Within KeepPairs \
        --outFilterScoreMinOverLread 0.3 \
        --outFilterMatchNminOverLread 0.3 \
        {custom_params} \
        --runThreadN {num_threads}"
    elif not fq1 and not fq2 and bam:
        cmd = f"STAR --outFileNamePrefix {out_dir}/ \
        --runMode alignReads \
        --limitOutSAMoneReadBytes 1000000 \
        --genomeDir {index_dir} \
        --readFilesType SAM PE \
        --readFilesCommand 'samtools view' \
        --readFilesIn {bam} \
        --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
        --outSAMmode Full \
        --alignEndsType EndToEnd \
        --outFilterMultimapNmax -1 \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMunmapped Within KeepPairs \
        --outFilterScoreMinOverLread 0.3 \
        --outFilterMatchNminOverLread 0.3 \
        {custom_params} \
        --runThreadN {num_threads}"
    return cmd


def run(fq1, fq2, bam, index_dir, out_dir, threads, method, params):
    """Build up and run alignment command for the specified aligner.

    Args:
        fq1 (str): Input FQ1 file
        fq2 (str): Input FQ2 file
        bam (str): Input BAM file
        index_dir (str): Index path to use for processing
        out_dir (str): Output directory where results are stored to
        threads (int): Number of cpus to use for processing
        params (str): Custom parameters string to add to the final command

    Returns:
        str: command to be executed in the subprocess
    """
    cmd = None
    if method == "bowtie2":
        cmd = get_align_cmd_bowtie2(fq1, fq2, bam, index_dir, out_dir, threads, params)
    elif method == "star":
        cmd = get_align_cmd_star(fq1, fq2, bam, index_dir, out_dir, threads, params)
    #print(cmd)
    IOMethods.execute_cmd(cmd)
    #subprocess.run(cmd, shell=True, cwd=out_dir, check=None)

def add_aligner_args(parser):
    """Add parser arguments.

    Args:
        parser (obj): Command line parser object to use
    """
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
        "--output_dir",
        dest="output_dir",
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


def alignment_command(args) -> None:
    """Run command with command line arguments.

    Args:
        args (obj): Command line arguments
    """
    run(
        fq1=args.fq1,
        fq2=args.fq2,
        bam=args.bam,
        index_dir=args.index_dir,
        out_dir=args.output_dir,
        threads=args.threads,
        method=args.method,
        params=args.params
    )
