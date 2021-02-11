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


def get_fastq_files(paths):
    """This function returns a list of fastq files for the given list of paths."""
    fastqs = []
    for path in paths:
        if os.path.isdir(path):
            files = os.listdir(path)
            for filename in files:
                file_path = os.path.join(path, filename)
                if os.path.isfile(file_path) and filename.endswith((".fastq.gz", ".fq.gz", ".fq", ".fastq")):
                    fastqs.append(file_path)
                elif os.path.isdir(file_path):
                    fastqs_tmp = get_fastq_files([file_path])
                    fastqs.extend(fastqs_tmp)
        elif os.path.isfile(path) and path.endswith((".fastq.gz", ".fq.gz", ".fq", ".fastq")):
            fastqs.append(path)
    return fastqs


def pair_fastq_files(fastqs):
    """Pairs paired-end fastq files, if necessary using regular expressions."""
    left = []
    right = []
    sample_id = []
    # iterate over the sorted list of file names and check for left/right pair file

    #print("\nGoing to process the following read files...")
    for i, fq_file in enumerate(sorted(fastqs)):
        try:
            # Search for 1 or 2 between "_R|_" and "_|.f" in the basename of the file
            forrev = re.search('[_|.]R([1-2])(.*)(_|[.]f)', os.path.basename(fq_file)).group(1)
            #print(forrev)
        except AttributeError:
            forrev = '-1'

        if forrev == '1':
            left.append(fq_file)
        elif forrev == '2':
            right.append(fq_file)
        else:
            print('Warning: Ignoring \"{}\" as it doesn\'t seem to be a valid fastq file'.format(fq_file))
    # Check whether file names match between paired files
    if right:
        for i, _ in enumerate(left):
            sample_id_l = re.search('(.*)[_|.]R1(.*)(_|[.]f)', os.path.basename(left[i])).group(1)
            sample_id_r = re.search('(.*)[_|.]R2(.*)(_|[.]f)', os.path.basename(right[i])).group(1)
            sample_id.append(sample_id_l)

            if sample_id_l != sample_id_r:
                print('Error 99: Paired files names {0} and {1} do not match!'.format(sample_id_l, sample_id_r))
                sys.exit(99)
    return (left, right, sample_id)


class CustomTranscriptome(object):
    
    def __init__(self, cfg, input_paths, seq_tab, bp_distance, working_dir):
        
        self.cfg = cfg
        self.input_paths = []
        for path in input_paths:
            self.input_paths.append(path.rstrip("/"))
        self.seq_tab = seq_tab
        self.bp_distance = bp_distance
        self.working_dir = working_dir.rstrip("/")
        if not os.path.exists(self.working_dir):
            print("INFO: Creating folder: {}".format(self.working_dir))
            os.makedirs(self.working_dir)
        self.fasta = os.path.join(self.working_dir, "context.fa")
        self.bed = os.path.join(self.working_dir, "context.bed")
        self.csv_to_fasta()
        self.generate_bed()
        self.build_index()


    def calc_fasta_size(self):
        """This function calculates the FASTA size in bp."""
        fasta_size = 0
        with open(self.fasta) as inf:
            for line in inf:
                if not line.rstrip().startswith(">"):
                    fasta_size += len(line.rstrip())
        print("INFO: FASTA Size: {}bp".format(fasta_size))
        return fasta_size


    def csv_to_fasta(self):
        """This function converts the target sequences TSV/CSV file to the FASTA format."""
        outf_context = open(self.fasta, "w")
        
        with open(self.seq_tab, "r", newline="\n") as csvfile:
            # Auto detect dialect of input file
            dialect = csv.Sniffer().sniff(csvfile.readline(), delimiters=";,\t")
            csvfile.seek(0)
            reader = csv.DictReader(csvfile, dialect=dialect)

            # Iterate over input file rows
            for row in reader:
                
                name = row["name"]
                sequence = row["sequence"]
                position = row["position"]
                
                outf_context.write(">{}\n".format(name))
                for i in range(0, len(sequence), 60):
                    outf_context.write("{}\n".format(sequence[i:i+60]))
                    
        outf_context.close()

    def generate_bed(self):
        """This function converts the target sequences TSV/CSV file to the BED format."""
        outf_context = open(self.bed, "w")
        
        with open(self.seq_tab) as csvfile:
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

                outf_context.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(name, start, end1, name_pos + "_A", "0", "+"))
                outf_context.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(name, end1, end2, name_pos + "_B", "0", "+"))

        outf_context.close()

    def build_index(self):
        """This function builds an index on the input sequences for mapping with STAR."""
        
        star_idx_path = os.path.join(self.working_dir, "STAR_idx")


        for folder in [star_idx_path]:
            if not os.path.exists(folder):
                print("INFO: Creating folder: {}".format(folder))
                os.makedirs(folder)            

        shell_script = os.path.join(self.working_dir, "generate_index.sh")
        out_shell = open(shell_script, "w")
        
        # calculate length of SA pre-indexing string.
        # according to STAR parameter --genomeSAindexNbases 
        # If genome size is to small, mapping is very slow, therfore use at least 4
        sa_index_nbases = min(14, max(4, int(math.log(self.calc_fasta_size()) / 2 - 1)))
        cmd_star = "{} --runMode genomeGenerate --runThreadN 12 --genomeSAindexNbases {} --genomeDir {} --genomeFastaFiles {}".format(self.cfg.get('commands','star_cmd'), sa_index_nbases, star_idx_path, self.fasta)

        out_shell.write("#!/bin/sh\n\n")
        out_shell.write("working_dir={}\n".format(self.working_dir))
        out_shell.write("echo \"Generating STAR index\"\n")
        out_shell.write("{}\n".format(cmd_star))
        out_shell.write("echo \"Processing done!\"\n")
        out_shell.close()

        if not os.path.exists(os.path.join(star_idx_path, "Genome")):
            self.execute_cmd(cmd_star)

    def run(self):
        """This function starts processing of the samples."""
        
        fastqs = get_fastq_files(self.input_paths)
        file_array = sorted(fastqs)
        print("INFO: Fastq files: {}".format(file_array))

        left, right, sample_id = pair_fastq_files(fastqs)

        for i, _ in enumerate(left):
            if len(left) == len(right):
                print("INFO: Processing Sample ID: {} (paired end)".format(sample_id[i]))
                print("INFO:    Sample 1: {}".format(left[i]))
                print("INFO:    Sample 2: {}".format(right[i]))
                self.execute_pipeline(left[i], right[i])


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
                --outFilterMultimapNmax -1 \
                --outSAMattributes NH HI AS nM NM MD \
                --outSAMunmapped None \
                --outFilterMismatchNoverLmax 0.05 \
                --outFilterMismatchNoverReadLmax 0.05 \
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

    
    ct = CustomTranscriptome(cfg, args.input_paths, args.seq_tab, args.bp_distance, args.output_folder)
    ct.run()
    
    script_call = "python {} {}".format(os.path.realpath(__file__), " ".join(sys.argv[1:]))
    
    outf = open(os.path.join(args.output_folder, "run_command.sh"), "w")
    outf.write("#!/bin/sh\n\n")
    outf.write(script_call)
    outf.close()

if __name__ == "__main__":
    main()
