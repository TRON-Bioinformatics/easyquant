#!/usr/bin/env python

import os
import re
import sys
import pwd,grp,stat
import time
import csv
import math
from argparse import ArgumentParser
from misc.config import Config
import subprocess

class CustomTranscriptome(object):
    
    def __init__(self, cfg, input, seq_tab, bp_distance, working_dir):
        
        self.cfg = cfg
        self.input = []
        for i in input:
            self.input.append(i.rstrip("/"))
        self.seq_tab = seq_tab
        self.bp_distance = bp_distance
        self.working_dir = working_dir.rstrip("/")
        if not os.path.exists(self.working_dir):
            print "INFO: Creating folder: " + self.working_dir
            os.makedirs(self.working_dir)
        self.fasta = os.path.join(self.working_dir,"context.fa")
        self.bed = os.path.join(self.working_dir,"context.bed")
        self.csv_to_fasta()
        self.generate_bed()
        self.build_index()


    def get_fastq_files(self, paths):
        '''This function returns a list of fastq files for the given list of paths.'''
        fastqs = []
        for path in paths:
            if os.path.isdir(path):
                files = os.listdir(path)
                for file in files:
                    file_path = os.path.join(path,file)
                    if os.path.isfile(file_path) and file.endswith((".fastq.gz",".fq.gz",".fq",".fastq")):
                        fastqs.append(file_path)
                    elif os.path.isdir(file_path):
                        fastqs_tmp = self.get_fastq_files([file_path])
                        fastqs.extend(fastqs_tmp)
            elif os.path.isfile(path) and path.endswith((".fastq.gz",".fq.gz",".fq",".fastq")):
                fastqs.append(path)
        return fastqs

    def calc_genome_size(self):
        genome_size = 0
        with open(self.fasta) as inf:
            for line in inf:
                if not line.rstrip().startswith(">"):
                    genome_size += len(line.rstrip())
        print "INFO: Custom Genome Size: " + str(genome_size) + "bp"
        return genome_size

### FGID;Fusion_Gene;Breakpoint1;Breakpoint2;FTID;context_sequence_id;type;exon_boundary1;exon_boundary2;exon_boundary;bp1_frame;bp2_frame;frame;wt1_expression;wt2_expression;context_sequence;context_sequence_bp;wt1_context_sequence;wt1_context_sequence_bp;wt2_context_sequence;wt2_context_sequence_bp;neo_peptide_sequence;neo_peptide_sequence_bp;transcript_sequence ###
### 
    def csv_to_fasta(self):
        
        outf_context = open(self.fasta,"w")
        
        with open(self.seq_tab) as csvfile:
            inf = csv.reader(csvfile, delimiter=';')
            header = inf.next()

            col = {}
            for i,colname in enumerate(header):
                col[colname] = i
                
            for c,row in enumerate(inf):
                
                name = row[col["name"]]
                sequence = row[col["sequence"]]
                position = row[col["position"]]
                
                outf_context.write(">" + name + "\n")
                for i in range(0, len(sequence), 60):
                    outf_context.write(sequence[i:i+60]+"\n")
                    
        outf_context.close()

    def generate_bed(self):
        
        outf_context = open(self.bed,"w")
        
        with open(self.seq_tab) as csvfile:
            
            inf = csv.reader(csvfile, delimiter=';')
            header = inf.next()
            col = {}
            
            # read column names as dictionary to column numbers
            for i,colname in enumerate(header):
                col[colname] = i
            
            # iterate over rows
            for c,row in enumerate(inf):
                
                name = row[col["name"]]
                sequence = row[col["sequence"]]
                position = row[col["position"]]
                
                name_pos = name + "_pos_" + position
                
                # use 0 as position, if not given
                location = 0
                if position != "NA":
                    location = int(position)

                strand1 = "+"
                strand2 = "+"

                start = 0
                end1 = location
                end2 = len(sequence)

                outf_context.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (name, start, end1, name_pos + "_A", "0", "+"))
                outf_context.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (name, end1, end2,  name_pos + "_B", "0", "+"))

        outf_context.close()

    def build_index(self):
        '''This function builds an index on the input sequences for mapping with STAR.'''
        
        star_idx_path = os.path.join(self.working_dir,"STAR_idx")


        for folder in [star_idx_path]:
            if not os.path.exists(folder):
                print "INFO: Creating folder: " + folder
                os.makedirs(folder)            

        shell_script = os.path.join(self.working_dir,"generate_index.sh")
        out_shell = open(shell_script,"w")
        
        # calculate length of SA pre-indexing string.
        # according to STAR parameter --genomeSAindexNbases 
        # If genome size is to small, mapping is very slow, therfore use at least 4
        sa_index_nbases = min(14, max(4, int(math.log(self.calc_genome_size()) / 2 - 1)))
        cmd_star = "%s --runMode genomeGenerate --runThreadN 12 --genomeSAindexNbases %s --genomeDir %s --genomeFastaFiles %s" % (self.cfg.get('commands','star_cmd'),sa_index_nbases,star_idx_path,self.fasta)

        out_shell.write("#!/bin/sh\n\n")
        out_shell.write("working_dir="+self.working_dir+"\n")
        out_shell.write("echo \"Generating STAR index\"\n")
        out_shell.write(cmd_star+"\n")
        out_shell.write("echo \"Processing done!\"\n")
        out_shell.close()

        if not os.path.exists(os.path.join(star_idx_path,"Genome")):
            self.execute_cmd(cmd_star)

    def run(self):
        '''This function starts processing of the samples.'''
        
        fastqs = self.get_fastq_files(self.input)
        file_array = sorted(fastqs)
        
        print "INFO: Fastq files:", file_array
        
        left = []
        right = []
        for i in range(len(file_array)):
            file_array_split = file_array[i].rsplit("_",2)
            skewer_split = ""
            try:
                skewer_split = file_array[i].rsplit("-",1)[1].rsplit(".",2)[0]
            except:
                skewer_split = ""
            if "R1" in file_array_split[1] or skewer_split == "pair1":
                left.append(file_array[i])
            elif "R2" in file_array_split[1] or skewer_split == "pair2":
                right.append(file_array[i])
                
        for i,ele in enumerate(left):

            if len(right) != 0:
                print "INFO: Start pipeline with:"
                print "INFO:    Sample 1: " + left[i]
                print "INFO:    Sample 2: " + right[i]
                self.execute_pipeline(left[i],right[i])

    def execute_pipeline(self,file_1,file_2):
        
        print "INFO: Execute pipeline on", file_1, file_2
        
        star_genome_path = os.path.join(self.working_dir,"STAR_idx")

        star_path = os.path.join(self.working_dir,"star")
        
        # create folder
        if not os.path.exists(star_path):
            print "INFO: Creating folder: " + star_path
            os.makedirs(star_path)
        
        # define bash script in working directory    
        shell_script = os.path.join(self.working_dir, "requant.sh")
        
        # start to write shell script to execute mapping cmd
        with open(shell_script, "w") as out_shell:

            bam_file = os.path.join(star_path,"Aligned.sortedByCoord.out.bam")
    
            cmd_star = "%s --outFileNamePrefix %s --limitOutSAMoneReadBytes 1000000 --genomeDir %s --readFilesCommand 'gzip -d -c -f' --readFilesIn %s %s --outSAMmode Full --outFilterMultimapNmax -1 --outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverReadLmax 0.02 --outSAMtype BAM SortedByCoordinate --runThreadN 12 && %s index %s" % (self.cfg.get('commands','star_cmd'), star_path + "/", star_genome_path, file_1, file_2, self.cfg.get('commands', 'samtools_cmd'), bam_file)
    
            cmd_class = "%s -i %s -b %s -d %s -o %s" % (self.cfg.get('commands', 'classification_cmd'), bam_file, self.bed, self.bp_distance, os.path.join(self.working_dir, "quantification.csv"))
    
            out_shell.write("#!/bin/sh\n\n")
            out_shell.write("fq1="+file_1+"\n")
            if file_2:
                out_shell.write("fq2="+file_2+"\n")
            out_shell.write("working_dir="+self.working_dir+"\n")
            out_shell.write("echo \"Starting pipeline...\"\n")
            if file_2:
                out_shell.write("echo \"Starting STAR\"\n")
                out_shell.write(cmd_star+"\n")
                out_shell.write("echo \"Starting quantification\"\n")
                out_shell.write(cmd_class+"\n")
            out_shell.write("echo \"Processing done!\"\n")

        if not os.path.exists(os.path.join(star_path, "Log.final.out")):
            self.execute_cmd(cmd_star)

        if not os.path.exists(os.path.join(self.working_dir, "quantification.csv")):
            self.execute_cmd(cmd_class)

        print "INFO: Processing complete for " + self.working_dir

    def execute_cmd(self, cmd, working_dir = "."):
        
        print(cmd)
        p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = working_dir, shell = True)
        (stdoutdata, stderrdata) = p.communicate()
        r = p.returncode
        if r != 0:
            print("Error: Command \"{}\" returned non-zero exit status".format(cmd))
            print(stderrdata)
            sys.exit(1)

def main():
    parser = ArgumentParser(description='Processing of demultiplexed FASTQs')

    parser.add_argument('-i', '--input', dest='input', nargs='+', help='Specify the fastq folder(s) or fastq file(s) to process.',required=True)
    parser.add_argument('-s', '--sequence_tab', dest='seq_tab', help='Specify the reference sequences as table with colums name, sequence, and position',required=True)
    parser.add_argument('-d', '--bp_distance', dest='bp_distance', help='Threshold in base pairs for the required overlap size of reads on both sides of the breakpoint for junction/spanning read counting', default=10)
    parser.add_argument('-o', '--output-folder', dest='output_folder', help='Specify the folder to save the results into.',required=True)
    args = parser.parse_args()

    cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)),'config.ini'))

    
    ct = CustomTranscriptome(cfg, args.input, args.seq_tab, args.bp_distance, args.output_folder)
    ct.run()
    
    script_call = "python " + os.path.realpath(__file__) + " " + " ".join(sys.argv[1:])
    
    outf = open(os.path.join(args.output_folder,"run_command.sh"),"w")
    outf.write("#!/bin/sh\n\n")
    outf.write(script_call)
    outf.close()

if __name__ == '__main__':
    main()
