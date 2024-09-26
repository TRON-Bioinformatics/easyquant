"""
This module generates the indexing commands for the pipeline.
"""

import math
import subprocess

# pylint: disable=E0401
from Bio import SeqIO # type: ignore

def get_sequence_count_and_len(fasta_file) -> tuple:
    """Returns the number of sequences and the total bases from a fasta file."""
    seq_num = 0
    total_bases = 0
    with open(fasta_file, encoding="utf8") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            seq_num += 1
            total_bases += len(record.seq)
    return (seq_num, total_bases)



def get_index_cmd_bowtie2(fasta_in, index_out, num_threads) -> str:
    """_summary_

    Args:
        fasta_in (_type_): _description_
        index_out (_type_): _description_
        num_threads (_type_): _description_

    Returns:
        str: _description_
    """
    cmd = f"bowtie2-build --threads {num_threads} {fasta_in} {index_out}/bowtie"
    return cmd


def get_index_cmd_star(fasta_in, index_out, num_threads) -> str:
    """_summary_

    Args:
        fasta_in (_type_): _description_
        index_out (_type_): _description_
        num_threads (_type_): _description_

    Returns:
        str: _description_
    """
    (seq_num, fasta_size) = get_sequence_count_and_len(fasta_in)
    chr_bin_n_bits = min(18, int(math.log(fasta_size / seq_num, 2)))
    sa_index_nbases = min(14, max(4, int(math.log(fasta_size) / 2 - 1)))
    cmd = f"STAR --runMode genomeGenerate \
    --limitGenomeGenerateRAM 40000000000 \
    --runThreadN {num_threads} \
    --genomeChrBinNbits {chr_bin_n_bits} \
    --genomeSAindexNbases {sa_index_nbases} \
    --genomeDir {index_out} \
    --genomeFastaFiles {fasta_in}"

    return cmd


def run(fasta_in, index_out, threads, method) -> None:
    """_summary_

    Args:
        fasta_in (_type_): _description_
        index_out (_type_): _description_
        threads (_type_): _description_
        method (_type_): _description_
    """
    if method == "bowtie2":
        subprocess.run(get_index_cmd_bowtie2(fasta_in, index_out, threads).split(" "), check=False)
    elif method == "star":
        subprocess.run(get_index_cmd_star(fasta_in, index_out, threads).split(" "), check=False)


def add_indexing_args(parser):
    """_summary_

    Args:
        parser (_type_): _description_
    """
    parser.add_argument(
        "-i",
        "--input_fasta",
        dest="input_fasta",
        help="Specify input FASTA file",
        required=True
    )
    parser.add_argument(
        "-o",
        "--index_dir",
        dest="index_dir",
        help="Specify path to output index dir",
        default="custom_index"
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
    parser.set_defaults(func=indexing_command)


def indexing_command(args):
    """_summary_

    Args:
        args (_type_): _description_
    """
    run(args.input_fasta, args.index_dir, args.threads, args.method)
