#!/usr/bin/env python

from argparse import ArgumentParser
import csv
import grp
import math
import os
import pwd
import re
import stat
import subprocess
import sys
import time


from misc.config import Config
import misc.io_methods as IOMethods
from custom_transcriptome import CustomTranscriptome



class Easyquant(object):
    
    def __init__(self, cfg, input_paths, seq_tab, bp_distance, working_dir, interval_mode):
        
        self.cfg = cfg
        self.input_paths = []
        for path in input_paths:
            self.input_paths.append(path.rstrip("/"))
        self.seq_tab = seq_tab
        self.bp_distance = bp_distance
        self.working_dir = working_dir.rstrip("/")
        IOMethods.create_folder(self.working_dir)
        self.interval_mode = interval_mode


    def run(self, method):
        """This function starts processing of the samples."""

        ct = CustomTranscriptome(self.seq_tab, self.working_dir)
        ct.run(method)
        
        fastqs = IOMethods.get_fastq_files(self.input_paths)
        file_array = sorted(fastqs)
        print("INFO: Fastq files: {}".format(file_array))

        sample_list = IOMethods.pair_fastq_files(fastqs)

        for sample in sample_list:
            print("INFO: Processing Sample ID: {} (paired end)".format(sample[0]))
            print("INFO:    Sample 1: {}".format(sample[1]))
            print("INFO:    Sample 2: {}".format(sample[2]))
            self.execute_pipeline(sample[1], sample[2], method)


    def execute_pipeline(self, file_1, file_2, method):
        """This function runs the pipeline on paired-end FASTQ files."""
        print("INFO: Execute pipeline on {} {}".format(file_1, file_2))
        


        output_path = os.path.join(self.working_dir, method)
        
        # create folder

        if not os.path.exists(output_path):
            print("INFO: Creating folder: {}".format(output_path))
            os.makedirs(output_path)

        
        # define bash script in working directory    
        shell_script = os.path.join(self.working_dir, "requant.sh")
        cmd = ""
        # start to write shell script to execute mapping cmd
        with open(shell_script, "w") as out_shell:
            sam_file = os.path.join(output_path, "Aligned.out.sam")
            bam_file = os.path.join(output_path, "Aligned.sortedByCoord.out.bam")
            quant_file = os.path.join(self.working_dir, "quantification.tsv")
            #read_file = os.path.join(self.working_dir, "read_info.tsv")

            if method == "bowtie2":
                genome_path = os.path.join(self.working_dir, "{}_idx".format(method))
                cmd = "{0} -x {1}/bowtie -1 {2} -2 {3} -S {4}".format(
                    self.cfg.get('commands', 'bowtie2_cmd'),
                    genome_path,
                    file_1,
                    file_2,
                    sam_file
                )

            elif method == "bwa":
                genome_path = os.path.join(self.working_dir, "context.fa")
                cmd = "{0} mem {1} {2} {3} > {4}".format(
                    self.cfg.get('commands', 'bwa_cmd'), 
                    genome_path, 
                    file_1, 
                    file_2, 
                    sam_file
                )

            elif method == "star":
                genome_path = os.path.join(self.working_dir, "{}_idx".format(method))

                cmd = "{} --outFileNamePrefix {} \
                --limitOutSAMoneReadBytes 1000000 \
                --genomeDir {} \
                --readFilesCommand 'gzip -d -c -f' \
                --readFilesIn {} {} \
                --outSAMmode Full \
                --alignEndsType EndToEnd \
                --outFilterMultimapNmax -1 \
                --outSAMattributes NH HI AS nM NM MD \
                --outSAMunmapped None \
                --outFilterMismatchNoverLmax 0.05 \
                --outFilterMismatchNoverReadLmax 0.05 \
                --runThreadN 12".format(self.cfg.get('commands','star_cmd'), 
                                        output_path + "/", genome_path, 
                                        file_1, file_2)



            if self.interval_mode:
                cmd_class = "{} -i {} -t {} -d {} -o {} --interval-mode".format(
                    self.cfg.get('commands', 'classification_cmd'),
                    sam_file,
                    self.seq_tab,
                    self.bp_distance,
                    self.working_dir
                )
            else:
                cmd_class = "{} -i {} -t {} -d {} -o {}".format(
                    self.cfg.get('commands', 'classification_cmd'),
                    sam_file,
                    self.seq_tab,
                    self.bp_distance,
                    self.working_dir
                )
    
            out_shell.write("#!/bin/sh\n\n")
            out_shell.write("fq1={}\n".format(file_1))
            if file_2:
                out_shell.write("fq2={}\n".format(file_2))
            out_shell.write("working_dir={}\n".format(self.working_dir))
            out_shell.write("echo \"Starting pipeline...\"\n")
            if file_2:
                out_shell.write("echo \"Starting alignment\"\n")
                out_shell.write("{}\n".format(cmd))
                out_shell.write("echo \"Starting quantification\"\n")
                out_shell.write("{}\n".format(cmd_class))
            out_shell.write("echo \"Processing done!\"\n")

        if not os.path.exists(sam_file):
            self.execute_cmd(cmd)

        if not os.path.exists(quant_file):
            self.execute_cmd(cmd_class)

        print("INFO: Processing complete for {}".format(self.working_dir))

    def execute_cmd(self, cmd, working_dir = "."):
        """This function pushes a command into a subprocess."""
        print(cmd)
        p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = working_dir, shell = True)
        (stdoutdata, stderrdata) = p.communicate()
        r = p.returncode
        if r != 0:
            print("Error: Command \"{}\" returned non-zero exit status".format(cmd))
            print(stderrdata)
            sys.exit(1)

def main():
    parser = ArgumentParser(description="Processing of demultiplexed FASTQs")

    parser.add_argument("-i", "--input", dest="input_paths", nargs="+", help="Specify the fastq folder(s) or fastq file(s) to process.", required=True)
    parser.add_argument("-s", "--sequence_tab", dest="seq_tab", help="Specify the reference sequences as table with colums name, sequence, and position", required=True)
    parser.add_argument("-d", "--bp_distance", dest="bp_distance", help="Threshold in base pairs for the required overlap size of reads on both sides of the breakpoint for junction/spanning read counting", default=10)
    parser.add_argument("-o", "--output-folder", dest="output_folder", help="Specify the folder to save the results into.", required=True)
    parser.add_argument("-m", "--method", dest="method", choices=["star", "bowtie2", "bwa"], help="Specify alignment software to generate the index", default="star")
    parser.add_argument('--interval-mode', dest='interval_mode', action='store_true', help='Specify if interval mode shall be used')
    args = parser.parse_args()

    cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini"))

    eq = Easyquant(cfg, args.input_paths, args.seq_tab, args.bp_distance, args.output_folder, args.interval_mode)
    eq.run(args.method)

    script_call = "python {} {}".format(os.path.realpath(__file__), " ".join(sys.argv[1:]))
    
    outf = open(os.path.join(args.output_folder, "run_command.sh"), "w")
    outf.write("#!/bin/sh\n\n")
    outf.write(script_call)
    outf.close()

if __name__ == "__main__":
    main()
