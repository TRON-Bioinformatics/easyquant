#!/usr/bin/env python

import os
import re
import sys
import pwd,grp,stat
import time
import csv
import math
from argparse import ArgumentParser
import misc.io_methods as IOMethods
from misc.config import Config
from misc.samples import Samples
import misc.queue as Queueing


class CustomTranscriptome(object):
    
    def __init__(self, cfg, input, fusions_table, working_dir, partitions):
        
        csv.field_size_limit(999999999)
        self.cfg = cfg
        self.input = []
        for i in input:
            self.input.append(i.rstrip("/"))
        self.fusions_table = fusions_table
        self.working_dir = working_dir.rstrip("/")
        IOMethods.create_folder(self.working_dir)
        self.samples = Samples(os.path.join(self.working_dir,"samples.csv"))
        self.fasta = os.path.join(self.working_dir,"context.fa")
        self.bed = os.path.join(self.working_dir,"context.bed")
        self.partitions = partitions
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
        print "Custom Genome Size: " + str(genome_size) + "bp"
        return genome_size

### FGID;Fusion_Gene;Breakpoint1;Breakpoint2;FTID;context_sequence_id;type;exon_boundary1;exon_boundary2;exon_boundary;bp1_frame;bp2_frame;frame;wt1_expression;wt2_expression;context_sequence;context_sequence_bp;wt1_context_sequence;wt1_context_sequence_bp;wt2_context_sequence;wt2_context_sequence_bp;neo_peptide_sequence;neo_peptide_sequence_bp;transcript_sequence ###
### 
    def csv_to_fasta(self):
        outf_context = open(self.fasta,"w")
        with open(self.fusions_table) as csvfile:
            inf = csv.reader(csvfile, delimiter=';')
            header = inf.next()

            col = {}
            for i,colname in enumerate(header):
                col[colname] = i
            for c,row in enumerate(inf):
                name = row[col["name"]]
                sequence = row[col["sequence"]]
                position = row[col["position"]]
                # context_id = row[col["context_sequence_id"]]
                # context_seq = row[col["context_sequence"]].upper()
                
                outf_context.write(">" + name + "\n")
                for i in range(0, len(sequence), 60):
                    outf_context.write(sequence[i:i+60]+"\n")
        outf_context.close()

    def generate_bed(self):
        
        outf_context = open(self.bed,"w")
        
        with open(self.fusions_table) as csvfile:
            
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
                
                name_pos = name + pos
                
                # genes = fgid.split("_", 1)
                # gene1 = genes[0]
                # gene2 = genes[1]
                # 
                # rel_bp = row[col["context_sequence_bp"]]
                
                # use 0 as position, if not given
                location = 0
                if position != "NA":
                    location = int(position)

                strand1 = "+"
                strand2 = "+"

                start = 0
                end1 = location
                end2 = len(sequence)

                outf_context.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, start, end1, name_pos + "_A", "0", "+"))
                outf_context.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, end1, end2,  name_pos + "_B", "0", "+"))

        outf_context.close()

    def build_index(self):
        star_idx_path = os.path.join(self.working_dir,"STAR_idx")

        for folder in [star_idx_path]:
            IOMethods.create_folder(folder)

        shell_script = os.path.join(self.working_dir,"generate_index.sh")
        out_shell = open(shell_script,"w")
        sa_index_nbases = int(math.log(self.calc_genome_size()) / 2 - 1)
        cmd_star = "%s --runMode genomeGenerate --runThreadN 12 --genomeSAindexNbases %s --genomeDir %s --genomeFastaFiles %s" % (self.cfg.get('commands','star_cmd'),sa_index_nbases,star_idx_path,self.fasta)

        out_shell.write("#!/bin/sh\n\n")
        out_shell.write("working_dir="+self.working_dir+"\n")
        out_shell.write("echo \"Generating STAR index\"\n")
        out_shell.write(cmd_star+"\n")
        out_shell.write("echo \"Processing done!\"\n")
        out_shell.close()

        if not os.path.exists(os.path.join(star_idx_path,"Genome")):
            self.submit_job_nonqueue(cmd_star,star_idx_path)

    def run(self):
        '''This function starts processing of the samples.'''
        fastqs = self.get_fastq_files(self.input)
        file_array = sorted(fastqs)
#        print file_array
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
                print "Sample 1: " + left[i]
                print "Sample 2: " + right[i]
                self.execute_pipeline(left[i],right[i])

    def execute_pipeline(self,file_1,file_2):

        star_genome_path = os.path.join(self.working_dir,"STAR_idx")

        star_path = os.path.join(self.working_dir,"star")

        for folder in [star_path]:
            IOMethods.create_folder(folder)
        shell_script = os.path.join(self.working_dir,"requant.sh")
        out_shell = open(shell_script,"w")

        bam_file = os.path.join(star_path,"Aligned.sortedByCoord.out.bam")

        cmd_star = "%s --limitOutSAMoneReadBytes 1000000 --genomeDir %s --readFilesCommand 'gzip -d -c -f' --readFilesIn %s %s --outSAMmode Full --outFilterMultimapNmax -1 --outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverLmax 0.02 --outSAMtype BAM SortedByCoordinate --runThreadN 12 && %s index %s" % (self.cfg.get('commands','star_cmd'), star_genome_path, file_1, file_2, self.cfg.get('commands', 'samtools_cmd'), bam_file)

        cmd_class = "%s -i %s -b %s -o %s" % (self.cfg.get('commands', 'classification_cmd'), bam_file, os.path.join(self.working_dir, self.bed), os.path.join(self.working_dir, "Classification.csv"))

        out_shell.write("#!/bin/sh\n\n")
        out_shell.write("fq1="+file_1+"\n")
        if file_2:
            out_shell.write("fq2="+file_2+"\n")
        out_shell.write("working_dir="+self.working_dir+"\n")
        out_shell.write("echo \"Starting pipeline...\"\n")
        if file_2:
            out_shell.write("echo \"Starting STAR\"\n")
            out_shell.write(cmd_star+"\n")
            out_shell.write("echo \"Starting Classification\"\n")
            out_shell.write(cmd_class+"\n")
        out_shell.write("echo \"Processing done!\"\n")
        out_shell.close()

        if not os.path.exists(os.path.join(star_path, "Log.final.out")):
            self.submit_job_nonqueue(cmd_star, star_path)

        if not os.path.exists(os.path.join(self.working_dir, "Classification.csv")):
            self.submit_job_nonqueue(cmd_class, self.working_dir)

        print "Processing complete for " + self.working_dir

    def submit_job_nonqueue(self,cmd,output_results_folder):
        Queueing._submit_nonqueue(cmd,output_results_folder)

    def submit_job(self,module,sample_id,cmd,cores,mem_usage,output_results_folder,dependencies):
        '''This function submits a job with the corresponding resources to slurm.'''
        q = Queue()
        job_name = "-".join([module,sample_id])
        already_running = self.get_jobs_with_name(job_name)
        if len(already_running) == 0:

            q.submit_slurm(job_name,cmd,cores,mem_usage,output_results_folder,dependencies,self.partitions)
            time.sleep(3)
        else:
            print job_name + " already running!"

        return job_name
       

def main():
    parser = ArgumentParser(description='Processing of demultiplexed FASTQs')

    parser.add_argument('-i', '--input', dest='input', nargs='+', help='Specify the fastq folder(s) or fastq file(s) to process.',required=True)
    parser.add_argument('-f', '--fusions', dest='fusions', help='Specify the reference fusions table to base index on.',required=True)
    parser.add_argument('-o', '--output-folder', dest='output_folder', help='Specify the folder to save the results into.',required=True)
    args = parser.parse_args()

    cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)),'config.ini'))

    
    ct = CustomTranscriptome(cfg, args.input, args.fusions, args.output_folder, args.partitions)
    ct.run()
    
    script_call = "python " + os.path.realpath(__file__) + " " + " ".join(sys.argv[1:])
    
    outf = open(os.path.join(args.output_folder,"custom_transcriptome.sh"),"w")
    outf.write("#!/bin/sh\n\n")
    outf.write(script_call)
    outf.close()

if __name__ == '__main__':
    main()
