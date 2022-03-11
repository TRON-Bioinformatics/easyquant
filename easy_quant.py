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
    
    def __init__(self, cfg, input_paths, seq_tab, bp_distance, working_dir):
        
        self.cfg = cfg
        self.input_paths = []
        for path in input_paths:
            self.input_paths.append(path.rstrip("/"))
        self.seq_tab = seq_tab
        self.bp_distance = bp_distance
        self.working_dir = working_dir.rstrip("/")
        IOMethods.create_folder(self.working_dir)


    def run(self):
        """This function starts processing of the samples."""

        ct = CustomTranscriptome(self.seq_tab, self.working_dir)
        ct.run()
        
        fastqs = IOMethods.get_fastq_files(self.input_paths)
        file_array = sorted(fastqs)
        print("INFO: Fastq files: {}".format(file_array))

        sample_list = IOMethods.pair_fastq_files(fastqs)

        for sample in sample_list:
            print("INFO: Processing Sample ID: {} (paired end)".format(sample[0]))
            print("INFO:    Sample 1: {}".format(sample[1]))
            print("INFO:    Sample 2: {}".format(sample[2]))
            self.execute_pipeline(sample[1], sample[2])


    def execute_pipeline(self, file_1, file_2):
        """This function runs the pipeline on paired-end FASTQ files."""
        print("INFO: Execute pipeline on {} {}".format(file_1, file_2))
        
        star_genome_path = os.path.join(self.working_dir, "STAR_idx")

        star_path = os.path.join(self.working_dir, "star")
        
        # create folder

        if not os.path.exists(star_path):
            print("INFO: Creating folder: {}".format(star_path))
            os.makedirs(star_path)

        
        # define bash script in working directory    
        shell_script = os.path.join(self.working_dir, "requant.sh")
        
        # start to write shell script to execute mapping cmd
        with open(shell_script, "w") as out_shell:

            bam_file = os.path.join(star_path, "Aligned.sortedByCoord.out.bam")
    
            cmd_star = "{} --outFileNamePrefix {} \
                --limitOutSAMoneReadBytes 1000000 \
                --genomeDir {} \
                --readFilesCommand 'gzip -d -c -f' \
                --readFilesIn {} {} \
                --outSAMmode Full \
                --outFilterMultimapNmax 1000 \
                --alignIntronMin 0 \
                --outSAMattributes NH HI AS nM NM MD \
                --outSAMunmapped Within \
                --outSAMtype BAM SortedByCoordinate \
                --runThreadN 12 && \
                {} index {}".format(self.cfg.get('commands','star_cmd'), star_path + "/",
                                star_genome_path, file_1, file_2,
                                self.cfg.get('commands', 'samtools_cmd'), bam_file)

            cmd_class = "{} -i {} -t {} -d {} -o {}".format(self.cfg.get('commands', 'classification_cmd'),
                                                        bam_file,
                                                        self.seq_tab,
                                                        self.bp_distance,
                                                        os.path.join(self.working_dir, "quantification.tsv"))
    
            out_shell.write("#!/bin/sh\n\n")
            out_shell.write("fq1={}\n".format(file_1))
            if file_2:
                out_shell.write("fq2={}\n".format(file_2))
            out_shell.write("working_dir={}\n".format(self.working_dir))
            out_shell.write("echo \"Starting pipeline...\"\n")
            if file_2:
                out_shell.write("echo \"Starting STAR\"\n")
                out_shell.write("{}\n".format(cmd_star))
                out_shell.write("echo \"Starting quantification\"\n")
                out_shell.write("{}\n".format(cmd_class))
            out_shell.write("echo \"Processing done!\"\n")

        if not os.path.exists(os.path.join(star_path, "Log.final.out")):
            self.execute_cmd(cmd_star)

        if not os.path.exists(os.path.join(self.working_dir, "quantification.tsv")):
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
    args = parser.parse_args()

    cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini"))

    eq = Easyquant(cfg, args.input_paths, args.seq_tab, args.bp_distance, args.output_folder)
    eq.run()

    script_call = "python {} {}".format(os.path.realpath(__file__), " ".join(sys.argv[1:]))
    
    outf = open(os.path.join(args.output_folder, "run_command.sh"), "w")
    outf.write("#!/bin/sh\n\n")
    outf.write(script_call)
    outf.close()

if __name__ == "__main__":
    main()
