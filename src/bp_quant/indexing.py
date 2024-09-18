import math
import subprocess

from Bio import SeqIO

def get_sequence_count_and_len(fasta_file):
    """Returns the number of sequences and the total bases from a fasta file"""
    seq_num = 0
    total_bases = 0
    with open(fasta_file) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            seq_num += 1
            total_bases += len(record.seq)
    return seq_num, total_bases



def get_index_cmd_bowtie2(fasta_in, index_out, num_threads):
    cmd = "bowtie2-build --threads {0} {1} {2}/bowtie".format(
        num_threads,
        fasta_in,
        index_out
    )
    return cmd


def get_index_cmd_star(fasta_in, index_out, num_threads):
    (seq_num, fasta_size) = get_sequence_count_and_len(fasta_in)
    chr_bin_n_bits = min(18, int(math.log(fasta_size / seq_num, 2)))
    sa_index_nbases = min(14, max(4, int(math.log(fasta_size) / 2 - 1)))
    cmd = "STAR --runMode genomeGenerate \
    --limitGenomeGenerateRAM 40000000000 \
    --runThreadN {0} \
    --genomeChrBinNbits {1} \
    --genomeSAindexNbases {2} \
    --genomeDir {3} \
    --genomeFastaFiles {4}".format(
        num_threads,
        chr_bin_n_bits,
        sa_index_nbases,
        index_out,
        fasta_in
    )

    return cmd


def run(fasta_in, index_out, threads, method):
    if method == "bowtie2":
        subprocess.run(get_index_cmd_bowtie2(fasta_in, index_out, threads).split(" "))
    elif method == "star":
        subprocess.run(get_index_cmd_star(fasta_in, index_out, threads).split(" "))


def add_indexing_args(parser):
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
    run(args.input_fasta, args.index_dir, args.threads, args.method)
