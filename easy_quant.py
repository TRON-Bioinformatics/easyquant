#!/usr/bin/env python

from argparse import ArgumentParser
from configparser import ConfigParser
#import gzip
import logging
import math
import os
import subprocess
import sys

import io_methods as IOMethods


#def get_read_count(fq_file):
#    n = 0
#    with gzip.open(fq_file) as inf:
#        for line in inf:
#            n += 1
#    print(n)
#    return n/4

#def get_read_count(fq_file):
#    """Parses input FASTQ to get read count"""
#    ps = subprocess.Popen(("zcat", fq_file), stdout=subprocess.PIPE)
#    result = subprocess.check_output(("wc", "-l"), stdin=ps.stdout)
#    return int(result) / 4
    


class Easyquant(object):
    
    def __init__(self, input_files, input_format, seq_tab, bp_distance, working_dir, allow_mismatches, interval_mode):

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
        self.input_files = [os.path.abspath(x) for x in input_files]
        self.input_format = input_format
        self.seq_tab = os.path.abspath(seq_tab)
        self.bp_distance = bp_distance
        self.allow_mismatches = allow_mismatches
        self.interval_mode = interval_mode


    def run(self, method, num_threads, star_cmd_params):
        """This function runs the pipeline on paired-end FASTQ files."""


        script_call = "python {} {}".format(os.path.realpath(__file__), " ".join(sys.argv[1:]))
        
        with open(os.path.join(self.working_dir, "run_command.sh"), "w") as outf:
            outf.write("#!/bin/sh\n\n")
            outf.write("version=\"{}\"\n".format(self.cfg.get("general", "version")))
            outf.write("python {} \\\n".format(os.path.realpath(__file__)))
            outf.write("-i {} \\\n".format(" ".join(self.input_files)))
            outf.write("-f {} \\\n".format(self.input_format))
            outf.write("-s {} \\\n".format(self.seq_tab))
            outf.write("-o {} \\\n".format(self.working_dir))
            outf.write("-d {} \\\n".format(self.bp_distance))
            outf.write("-m {} \\\n".format(method))
            outf.write("-t {} \\\n".format(num_threads))
            if self.interval_mode:
                outf.write("--interval_mode \\\n")
            if self.allow_mismatches:
                outf.write("--allow_mismatches \\\n")
            if star_cmd_params:
                outf.write(star_cmd_params)


        logging.info("Executing easyquant {}".format(self.cfg.get("general", "version")))
        if self.input_format == "fastq":
            logging.info("FQ1={}".format(self.input_files[0]))
            logging.info("FQ2={}".format(self.input_files[1]))
        elif self.input_format == "bam":
            logging.info("BAM={}".format(self.input_files[0]))

        genome_path = os.path.join(self.working_dir, "index")
        align_path = os.path.join(self.working_dir, "alignment")

        fasta_file = os.path.join(self.working_dir, "context.fa")
        sam_file = os.path.join(align_path, "Aligned.out.sam")
        bam_file = ""
        if self.input_format == "fastq":
            bam_file = os.path.join(align_path, "Aligned.out.bam")
        elif self.input_format == "bam":
            bam_file = os.path.join(align_path, "Processed.out.bam")
        quant_file = os.path.join(self.working_dir, "quantification.tsv")
        #num_reads_file = os.path.join(self.working_dir, "Star_org_input_reads.txt")
        
        #create folders
        IOMethods.create_folder(genome_path)
        IOMethods.create_folder(align_path)

        IOMethods.csv_to_fasta(self.seq_tab, fasta_file)


        #if not os.path.exists(num_reads_file) or os.stat(num_reads_file).st_size == 0:
        #    with open(num_reads_file, "w") as outf:
        #        outf.write(str(int(get_read_count(self.fq1))))


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
                self.input_files[0],
                self.input_files[1],
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
                self.input_files[0], 
                self.input_files[1], 
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

            # Most useful STAR parameters and their default values:
            #--outFilterMismatchNoverReadLmax 0.3 \
            #--scoreDelOpen -2 \
            #--scoreInsOpen -2 \
            #--scoreDelBase -2 \
            #--scoreInsBase -2 \


            if self.input_format == "fastq":

                align_cmd = "{0} --outFileNamePrefix {1} \
                --limitOutSAMoneReadBytes 1000000 \
                --genomeDir {2} \
                --readFilesCommand 'gzip -d -c -f' \
                --readFilesIn {3} \
                --outSAMmode Full \
                --alignEndsType EndToEnd \
                --outFilterMultimapNmax -1 \
                --outSAMattributes NH HI AS nM NM MD \
                --outSAMunmapped None \
                {4} \
                --runThreadN {5}".format(
                    self.cfg.get('commands','star'), 
                    align_path + "/", 
                    genome_path, 
                    " ".join(self.input_files),
                    star_cmd_params,
                    num_threads
                )

            elif self.input_format == "bam":
                align_cmd = "{0} --outFileNamePrefix {1} \
                --runMode alignReads \
                --limitOutSAMoneReadBytes 1000000 \
                --genomeDir {2} \
                --readFilesType SAM PE \
                --readFilesCommand 'samtools view' \
                --readFilesIn {3} \
                --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
                --outSAMmode Full \
                --alignEndsType EndToEnd \
                --outFilterMultimapNmax -1 \
                --outSAMattributes NH HI AS nM NM MD \
                --outSAMunmapped None \
                {4} \
                --runThreadN {5}".format(
                    self.cfg.get('commands','star'), 
                    align_path + "/", 
                    genome_path, 
                    " ".join(self.input_files),
                    star_cmd_params,
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

        #quant_cmd = "{0} -i {1} -d {2} -o {3}".format(
        #    self.cfg.get('commands', 'quantification'),
        #    bam_file,
        #    self.bp_distance,
        #    quant_file
        #)


        # define bash script in working directory    
        shell_script = os.path.join(self.working_dir, "requant.sh")
        # start to write shell script to execute mapping cmd
        with open(shell_script, "w") as out_shell:
            out_shell.write("#!/bin/sh\n\n")
            if self.input_format == "fastq":
                out_shell.write("fq1={}\n".format(self.input_files[0]))
                out_shell.write("fq2={}\n".format(self.input_files[1]))
            elif self.input_format == "bam":
                out_shell.write("bam={}\n".format(self.input_files[0]))
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

        if not os.path.exists(sam_file) or os.stat(sam_file).st_size == 0:
            IOMethods.execute_cmd(align_cmd)

        if not os.path.exists(bam_file) or os.stat(bam_file).st_size == 0:
            IOMethods.execute_cmd(sam_to_bam_cmd)

        if not os.path.exists(quant_file) or os.stat(quant_file).st_size == 0:
            IOMethods.execute_cmd(quant_cmd)


        logging.info("Processing complete for {}".format(self.working_dir))



def main():
    parser = ArgumentParser(description="Processing of demultiplexed FASTQs")

    parser.add_argument("-i", "--input_files", dest="input_files", nargs="+", help="Specify input file(s)", required=True)
    parser.add_argument("-f", "--input_format", dest="input_format", choices=["fastq", "bam"], help="Specify input format", default="fastq")
    parser.add_argument("-s", "--sequence_tab", dest="seq_tab", help="Specify the reference sequences as table with colums name, sequence, and position", required=True)
    parser.add_argument("-o", "--output-folder", dest="output_folder", help="Specify the folder to save the results into.", required=True)
    parser.add_argument("-d", "--bp_distance", dest="bp_distance", type=int, help="Threshold in base pairs for the required overlap size of reads on both sides of the breakpoint for junction/spanning read counting", default=10)
    parser.add_argument("--allow_mismatches", dest="allow_mismatches", action="store_true", help="Allow mismatches within the region around the breakpoint determined by the bp_distance parameter")
    parser.add_argument("--interval_mode", dest="interval_mode", action="store_true", help="Specify if interval mode shall be used")
    parser.add_argument("-m", "--method", dest="method", choices=["star", "bowtie2", "bwa"], help="Specify alignment software to generate the index", default="star")
    parser.add_argument("-t", "--threads", dest="num_threads", type=int, help="Specify number of threads to use for the alignment", default=1)
    parser.add_argument("--star_cmd_params", dest="star_cmd_params", help="Specify STAR commandline parameters to use for the alignment", default="")
    args = parser.parse_args()

    if args.input_format == "fastq" and len(args.input_files) != 2:
        print("FASTQ input format selected, but input requires two FASTQ files. Please try again.")
        sys.exit(1)

    if args.input_format == "bam" and len(args.input_files) > 1:
        print("BAM input format selected, but more than one file given. Please try again.")
        sys.exit(1)

    eq = Easyquant(args.input_files, args.input_format, args.seq_tab, args.bp_distance, args.output_folder, args.allow_mismatches, args.interval_mode)
    eq.run(args.method, args.num_threads, args.star_cmd_params)


if __name__ == "__main__":
    main()
