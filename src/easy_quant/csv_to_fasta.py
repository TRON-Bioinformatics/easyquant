import csv
import logging
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

csv.field_size_limit(sys.maxsize)

def csv_to_fasta(csv_in, fasta_out):
    """This function converts the target sequences TSV/CSV file to the FASTA format."""
    with open(fasta_out, "w") as fasta_handle, \
         open(csv_in, "r", newline="\n") as csvfile:
        # Auto detect dialect of input file
        dialect = csv.Sniffer().sniff(csvfile.readline(), delimiters=";,\t")
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, dialect=dialect)

        sequences = []
        # Iterate over input file rows
        for row in reader:
            record = SeqRecord(
                Seq(row["sequence"]),
                id=row["name"],
                name=row["name"],
                description=row["position"],
            )
            sequences.append(record)
        SeqIO.write(sequences, fasta_handle, "fasta")

    
def add_csv_to_fasta_args(parser):
    """Add csv_to_fasta arguments to parser"""
    parser.add_argument(
        "--input_csv", dest="in_csv", help="Specify input CSV file to convert", required=True
    )
    parser.add_argument(
        "--output_fasta", dest="out_fasta", help="Specify output FASTA file", required=True
    )
    parser.set_defaults(func=csv_to_fasta_command)


def csv_to_fasta_command(args):
    """Run csv_to_fasta conversion"""
    csv_to_fasta(args.in_csv, args.out_fasta)
