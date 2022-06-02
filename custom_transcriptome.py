#!/usr/bin/env python

from argparse import ArgumentParser

import csv
import math
import os
import subprocess
import sys

from misc.config import Config
import misc.io_methods as IOMethods

csv.field_size_limit(sys.maxsize)


def csv_to_fasta(csv_in, fasta_out):
    """This function converts the target sequences TSV/CSV file to the FASTA format."""
    outf = open(fasta_out, "w")


    with open(csv_in, "r", newline="\n") as csvfile:
        # Auto detect dialect of input file
        dialect = csv.Sniffer().sniff(csvfile.readline(), delimiters=";,\t")
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, dialect=dialect)

        # Iterate over input file rows
        for row in reader:

            name = row["name"]
            sequence = row["sequence"]
            position = row["position"]

            outf.write(">{}\n".format(name))
            for i in range(0, len(sequence), 60):
                outf.write("{}\n".format(sequence[i:i+60]))

    outf.close()


def get_fasta_size(fasta):
    """This function calculates the FASTA size in bp."""
    fasta_size = 0
    with open(fasta) as inf:
        for line in inf:
            if not line.rstrip().startswith(">"):
                fasta_size += len(line.rstrip())
    print("INFO: FASTA Size: {}bp".format(fasta_size))
    return fasta_size


def csv_to_bed(csv_in, bed_out):
    """This function converts the target sequences TSV/CSV file to the BED format."""
    outf = open(bed_out, "w")

    with open(csv_in) as csvfile:
        # Auto detect dialect of input file
        dialect = csv.Sniffer().sniff(csvfile.readline(), delimiters=";,\t")
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, dialect=dialect)

        # Iterate over input file rows
        for row in reader:

            name = row["name"]
            sequence = row["sequence"]
            position = row["position"]
            
            name_pos = "{}_pos_{}".format(name, position)

            # use 0 as position, if not given
            location = 0
            if position != "NA":
                location = int(position)

            strand1 = "+"
            strand2 = "+"

            start = 0
            end1 = location
            end2 = len(sequence)

            outf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(name, start, end1, name_pos + "_A", "0", "+"))
            outf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(name, end1, end2, name_pos + "_B", "0", "+"))

    outf.close()


def build_index(cfg, fasta_in, folder_out, method):
    """This function builds an index on the input sequences for mapping with BWA."""
    
    shell_script = os.path.join(folder_out, "generate_index.sh")
    out_shell = open(shell_script, "w")

    cmd = ""

    if method == "star":
        idx_path = os.path.join(folder_out, "star_idx")
        IOMethods.create_folder(idx_path)
        # calculate length of SA pre-indexing string.
        # according to STAR parameter --genomeSAindexNbases
        # If genome size is to small, mapping is very slow, therfore use at least 4
        sa_index_nbases = min(14, max(4, int(math.log(get_fasta_size(fasta_in)) / 2 - 1)))
        cmd = "{} --runMode genomeGenerate --limitGenomeGenerateRAM 40000000000 --runThreadN 12 --genomeSAindexNbases {} --genomeDir {} --genomeFastaFiles {}".format(cfg.get('commands','star_cmd'), sa_index_nbases, idx_path, fasta_in)
    elif method == "bwa":
        cmd = "{} index {}".format(cfg.get('commands', 'bwa_cmd'), fasta_in)
    elif method == "bowtie2":
        idx_path = os.path.join(folder_out, "bowtie2_idx")
        IOMethods.create_folder(idx_path)
        cmd = "{}-build {} {}/bowtie".format(cfg.get('commands', 'bowtie2_cmd'), fasta_in, idx_path)

    out_shell.write("#!/bin/sh\n\n")
    out_shell.write("working_dir={}\n".format(folder_out))
    out_shell.write("echo \"Generating {} index\"\n".format(method))
    out_shell.write("{}\n".format(cmd))
    out_shell.write("echo \"Processing done!\"\n")
    out_shell.close()

    #if not os.path.exists(os.path.join(star_idx_path, "Genome")):
    print("Executing CMD: {}".format(cmd))
    p = subprocess.run(cmd, shell = True)


class CustomTranscriptome(object):

    def __init__(self, seq_tab, working_dir):

        self.cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini"))
        self.working_dir = working_dir
        IOMethods.create_folder(self.working_dir)
        
        self.seq_tab = seq_tab
        self.fasta = os.path.join(self.working_dir, "context.fa")
        self.bed = os.path.join(self.working_dir, "context.bed")


    def run(self, method):
        csv_to_fasta(self.seq_tab, self.fasta)
        build_index(self.cfg, self.fasta, self.working_dir, method)


def main():

    parser = ArgumentParser(description="Generate STAR Index file for sequence CSV")
    parser.add_argument("-i", "--input_file", dest="input_file", help="Specifiy input CSV/TSV with sequence info")
    parser.add_argument("-o", "--output_folder", dest="output_folder", help="Specify output folder to store results into")
    parser.add_argument("-m", "--method", dest="method", choices=["star", "bowtie2", "bwa"], help="Specify alignment software to generate the index", default="star")
    args = parser.parse_args()

    seq_tab = args.input_file
    working_dir = args.output_folder

    ct = CustomTranscriptome(seq_tab, working_dir)
    ct.run(args.method)

if __name__ == "__main__":
    main()