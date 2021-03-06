#!/usr/bin/env python

from argparse import ArgumentParser
from configparser import ConfigParser
import logging
import math
import os
import subprocess
import sys

import io_methods as IOMethods



class Easyquant(object):
    
    def __init__(self, fq1, fq2, seq_tab, bp_distance, working_dir, allow_mismatches, interval_mode):

        self.working_dir = os.path.abspath(working_dir)
        IOMethods.create_folder(self.working_dir)

        logfile = os.path.join(self.working_dir, 'run.log')
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers=[
                logging.FileHandler(logfile),
                logging.StreamHandler()
            ]
        )


        cfg_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini")
        if not os.path.exists(cfg_file):
            logging.error("Config file does not exist! Please move and update your config file")
            logging.error("mv config.ini.sample config.ini -> edit")
            sys.exit(1)

        self.cfg = ConfigParser()
        self.cfg.read(cfg_file)
        self.fq1 = os.path.abspath(fq1)
        self.fq2 = os.path.abspath(fq2)
        self.seq_tab = os.path.abspath(seq_tab)
        self.bp_distance = bp_distance
        self.allow_mismatches = allow_mismatches
        self.interval_mode = interval_mode


    def run(self, method, num_threads):
        """This function runs the pipeline on paired-end FASTQ files."""


        script_call = "python {} {}".format(os.path.realpath(__file__), " ".join(sys.argv[1:]))
        
        with open(os.path.join(self.working_dir, "run_command.sh"), "w") as outf:
            outf.write("#!/bin/sh\n\n")
            outf.write("version=\"{}\"\n".format(self.cfg.get("general", "version")))
            outf.write("python {} \\\n".format(os.path.realpath(__file__)))
            outf.write("-1 {} \\\n".format(self.fq1))
            outf.write("-2 {} \\\n".format(self.fq2))
            outf.write("-s {} \\\n".format(self.seq_tab))
            outf.write("-o {} \\\n".format(self.working_dir))
            outf.write("-d {} \\\n".format(self.bp_distance))
            outf.write("-m {} \\\n".format(method))
            if self.interval_mode:
                outf.write("--interval-mode")


        logging.info("Executing easyquant {}".format(self.cfg.get("general", "version")))
        logging.info("FQ1={}".format(self.fq1))
        logging.info("FQ2={}".format(self.fq2))

        genome_path = os.path.join(self.working_dir, "index")
        align_path = os.path.join(self.working_dir, "alignment")

        fasta_file = os.path.join(self.working_dir, "context.fa")
        sam_file = os.path.join(align_path, "Aligned.out.sam")
        bam_file = os.path.join(align_path, "Aligned.out.bam")
        quant_file = os.path.join(self.working_dir, "quantification.tsv")

        
        #create folders
        IOMethods.create_folder(genome_path)
        IOMethods.create_folder(align_path)

        IOMethods.csv_to_fasta(self.seq_tab, fasta_file)

        index_cmd = None
        align_cmd = None
        quant_cmd = None


        if method == "bowtie2":
            index_cmd = "{0}-build {1} {2}/bowtie".format(
                self.cfg.get('commands', 'bowtie2'), 
                fasta_file, 
                genome_path
            )
            align_cmd = "{0} -p {1} -x {2}/bowtie -1 {3} -2 {4} -S {5}".format(
                self.cfg.get('commands', 'bowtie2'),
                num_threads,
                genome_path,
                self.fq1,
                self.fq2,
                sam_file
            )

        elif method == "bwa":
            genome_path = fasta_file
            index_cmd = "{0} index {1}".format(
                self.cfg.get('commands', 'bwa'), 
                fasta_file
            )
            align_cmd = "{0} mem -t {1} {2} {3} {4} > {5}".format(
                self.cfg.get('commands', 'bwa'), 
                num_threads,
                genome_path, 
                self.fq1, 
                self.fq2, 
                sam_file
            )

        elif method == "star":
            fasta_size = IOMethods.get_fasta_size(fasta_file)
            sa_index_nbases = min(14, max(4, int(math.log(fasta_size) / 2 - 1)))
            index_cmd = "{0} --runMode genomeGenerate \
            --limitGenomeGenerateRAM 40000000000 \
            --runThreadN {1} \
            --genomeSAindexNbases {2} \
            --genomeDir {3} \
            --genomeFastaFiles {4}".format(
                self.cfg.get('commands','star'), 
                num_threads,
                sa_index_nbases, 
                genome_path, 
                fasta_file
            )

            align_cmd = "{0} --outFileNamePrefix {1} \
            --limitOutSAMoneReadBytes 1000000 \
            --genomeDir {2} \
            --readFilesCommand 'gzip -d -c -f' \
            --readFilesIn {3} {4} \
            --outSAMmode Full \
            --alignEndsType EndToEnd \
            --outFilterMultimapNmax -1 \
            --outSAMattributes NH HI AS nM NM MD \
            --outSAMunmapped None \
            --outFilterMismatchNoverReadLmax {5} \
            --scoreDelOpen {6} \
            --scoreInsOpen {6} \
            --scoreDelBase {7} \
            --scoreInsBase {7} \
            --runThreadN {8}".format(
                self.cfg.get('commands','star'), 
                align_path + "/", 
                genome_path, 
                self.fq1, 
                self.fq2, 
                self.cfg.get('general', 'mismatch_ratio'),
                self.cfg.get('general', 'indel_open_penalty'),
                self.cfg.get('general', 'indel_extension_penalty'),
                num_threads
            )

        sam_to_bam_cmd = "{0} sort -o {2} {1} && {0} index {2}".format(
            self.cfg.get('commands', 'samtools'), 
            sam_file, 
            bam_file
        )

        allow_mismatches_str = ""
        interval_mode_str = ""
        if self.allow_mismatches:
            allow_mismatches_str = " --allow_mismatches"
        if self.interval_mode:
            interval_mode_str = " --interval_mode"

        quant_cmd = "{0} -i {1} -t {2} -d {3} -o {4}{5}{6}".format(
            self.cfg.get('commands', 'quantification'),
            sam_file,
            self.seq_tab,
            self.bp_distance,
            self.working_dir,
            allow_mismatches_str,
            interval_mode_str
        )


        # define bash script in working directory    
        shell_script = os.path.join(self.working_dir, "requant.sh")
        # start to write shell script to execute mapping cmd
        with open(shell_script, "w") as out_shell:
            out_shell.write("#!/bin/sh\n\n")
            out_shell.write("fq1={}\n".format(self.fq1))
            out_shell.write("fq2={}\n".format(self.fq2))
            out_shell.write("working_dir={}\n".format(self.working_dir))
            out_shell.write("echo \"Starting pipeline...\"\n")
            out_shell.write("echo \"Generating index\"\n")
            out_shell.write("{}\n".format(index_cmd))
            out_shell.write("echo \"Starting alignment\"\n")
            out_shell.write("{}\n".format(align_cmd))
            out_shell.write("echo \"Starting quantification\"\n")
            out_shell.write("{}\n".format(quant_cmd))
            out_shell.write("echo \"Processing done!\"\n")


        IOMethods.execute_cmd(index_cmd)

        if not os.path.exists(sam_file):
            IOMethods.execute_cmd(align_cmd)

        if not os.path.exists(quant_file):
            IOMethods.execute_cmd(quant_cmd)

        if not os.path.exists(bam_file):
            IOMethods.execute_cmd(sam_to_bam_cmd)

        logging.info("Processing complete for {}".format(self.working_dir))



def main():
    parser = ArgumentParser(description="Processing of demultiplexed FASTQs")

    parser.add_argument("-1", "--fq1", dest="fq1", help="Specify path to Read 1 (R1) FASTQ file", required=True)
    parser.add_argument("-2", "--fq2", dest="fq2", help="Specify path to Read 2 (R2) FASTQ file", required=True)
    parser.add_argument("-s", "--sequence_tab", dest="seq_tab", help="Specify the reference sequences as table with colums name, sequence, and position", required=True)
    parser.add_argument("-o", "--output-folder", dest="output_folder", help="Specify the folder to save the results into.", required=True)
    parser.add_argument("-d", "--bp_distance", dest="bp_distance", type=int, help="Threshold in base pairs for the required overlap size of reads on both sides of the breakpoint for junction/spanning read counting", default=10)
    parser.add_argument("--allow_mismatches", dest="allow_mismatches", action="store_true", help="Allow mismatches within the region around the breakpoint determined by the bp_distance parameter")
    parser.add_argument("--interval_mode", dest="interval_mode", action="store_true", help="Specify if interval mode shall be used")
    parser.add_argument("-m", "--method", dest="method", choices=["star", "bowtie2", "bwa"], help="Specify alignment software to generate the index", default="star")
    parser.add_argument("-t", "--threads", dest="num_threads", type=int, help="Specify number of threads to use for the alignment", default=1)
    args = parser.parse_args()

    eq = Easyquant(args.fq1, args.fq2, args.seq_tab, args.bp_distance, args.output_folder, args.allow_mismatches, args.interval_mode)
    eq.run(args.method, args.num_threads)


if __name__ == "__main__":
    main()
