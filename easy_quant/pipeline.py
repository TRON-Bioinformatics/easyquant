#!/usr/bin/env python

from argparse import ArgumentParser
from configparser import ConfigParser
import logging
import math
import os
import subprocess
import sys

import easy_quant.io_methods as IOMethods


def get_read_count(infile, tool, format="fq"):
    if format == "fq":
        return __get_read_count_fq(infile)
    elif format == "bam":
        return __get_read_count_bam(infile, tool)

def __get_read_count_fq(fq_file):
    """Parses input FASTQ to get read count"""
    ps = subprocess.Popen(("zcat", fq_file), stdout=subprocess.PIPE)
    result = subprocess.check_output(("wc", "-l"), stdin=ps.stdout)
    return int(result) / 2

def __get_read_count_bam(bam_file, tool):
    """Parses input BAM to get read count"""
    result = subprocess.check_output([tool, "view", "-c", bam_file])
    return int(result)


class Pipeline(object):
    
    def __init__(self, fq1, fq2, bam, seq_tab, bp_distance, working_dir, allow_mismatches, interval_mode, keep_bam, keep_all):

        self.working_dir = os.path.abspath(working_dir)
        self.module_dir = os.path.dirname(os.path.realpath(__file__))
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
        self.fq1 = None
        self.fq2 = None
        self.bam = None
        if fq1:
            self.fq1 = os.path.abspath(fq1)
        if fq2:
            self.fq2 = os.path.abspath(fq2)
        if bam:
            self.bam = os.path.abspath(bam)

        self.seq_tab = os.path.abspath(seq_tab)
        self.bp_distance = bp_distance
        self.allow_mismatches = allow_mismatches
        self.interval_mode = interval_mode
        self.keep_bam = keep_bam
        self.keep_all = keep_all

    def run(self, method, num_threads, star_cmd_params):
        """This function runs the pipeline on paired-end FASTQ files."""


        script_call = "python {} {}".format(os.path.realpath(__file__), " ".join(sys.argv[1:]))
        
        with open(os.path.join(self.working_dir, "run_command.sh"), "w") as outf:
            outf.write("#!/bin/sh\n\n")
            outf.write("version=\"{}\"\n".format(self.cfg.get("general", "version")))
            outf.write("python {} \\\n".format(os.path.realpath(__file__)))
            if self.fq1 and self.fq2:
                outf.write("-1 {} \\\n".format(self.fq1))
                outf.write("-2 {} \\\n".format(self.fq2))
            elif self.bam:
                outf.write("-b {} \\\n".format(self.bam))
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
        if self.fq1 and self.fq2:
            logging.info("FQ1={}".format(self.fq1))
            logging.info("FQ2={}".format(self.fq2))
        elif self.bam:
            logging.info("BAM={}".format(self.bam))

        genome_path = os.path.join(self.working_dir, "index")
        align_path = os.path.join(self.working_dir, "alignment")

        fasta_file = os.path.join(self.working_dir, "context.fa")
        sam_file = os.path.join(align_path, "Aligned.out.sam")
        bam_file = ""
        if self.fq1 and self.fq2:
            bam_file = os.path.join(align_path, "Aligned.out.bam")
        elif self.bam:
            bam_file = os.path.join(align_path, "Processed.out.bam")
        quant_file = os.path.join(self.working_dir, "quantification.tsv")
        num_reads_file = os.path.join(self.working_dir, "num_reads.txt")
        # Define files to be deleted after successful run
        clean_up_files = [genome_path, fasta_file, sam_file]
        if not self.keep_bam:
            clean_up_files.append(bam_file)
            clean_up_files.append("{}.bai".format(bam_file))
        
        #create folders
        IOMethods.create_folder(genome_path)
        IOMethods.create_folder(align_path)


        if not os.path.exists(num_reads_file) or os.stat(num_reads_file).st_size == 0:
            with open(num_reads_file, "w") as outf:
                if self.fq1:
                    outf.write(str(int(get_read_count(self.fq1, None, "fq"))))
                elif self.bam:
                    outf.write(str(int(get_read_count(self.bam, self.cfg.get('commands', 'samtools'), "bam"))))

        csv_to_fasta_cmd = "easy_quant csv2fasta \
        --input_csv {} \
        --output_fasta {}".format(
            self.seq_tab,
            fasta_file
        )

        index_cmd = "easy_quant index \
        --input_fasta {} \
        --index_dir {} \
        -t {} \
        --method {}".format(
            fasta_file,
            genome_path,
            num_threads,
            method
        )
        align_cmd = "easy_quant align --fq1 {} --fq2 {} --index_dir {} --output_path {} -t {} -m {}".format(
            self.fq1,
            self.fq2,
            genome_path,
            sam_file,
            num_threads,
            method
        )

        if method == "bwa":
            # Add index files from bwa
            clean_up_files.extend([
                "{}.amb".format(fasta_file),
                "{}.ann".format(fasta_file),
                "{}.bwt".format(fasta_file),
                "{}.pac".format(fasta_file),
                "{}.sa".format(fasta_file)
            ])

        elif method == "star":
            clean_up_files.extend([
                "{}/Log.progress.out".format(align_path),
                "{}/SJ.out.tab".format(align_path),
                "{}/_STARtmp".format(align_path),
            ])


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

        quant_cmd = "easy_quant count \
        -i {0} \
        -t {1} \
        -d {2} \
        -o {3}{4}{5}".format(
            sam_file,
            self.seq_tab,
            self.bp_distance,
            self.working_dir,
            allow_mismatches_str,
            interval_mode_str
        )

        clean_cmd = "for file in {}; \
            do {{ test -f $file && rm $file; }} || \
            {{ test -d $file && rm -r $file; }} done".format(" ".join(clean_up_files))

        # define bash script in working directory    
        shell_script = os.path.join(self.working_dir, "requant.sh")
        # start to write shell script to execute mapping cmd
        with open(shell_script, "w") as out_shell:
            out_shell.write("#!/bin/sh\n\n")
            if self.fq1 and self.fq2:
                out_shell.write("fq1={}\n".format(self.fq1))
                out_shell.write("fq2={}\n".format(self.fq2))
            elif self.bam:
                out_shell.write("bam={}\n".format(self.bam))
            out_shell.write("working_dir={}\n".format(self.working_dir))
            out_shell.write("echo \"Starting pipeline...\"\n")
            out_shell.write("echo \"Generating index\"\n")
            out_shell.write("{}\n".format(index_cmd))
            out_shell.write("echo \"Starting alignment\"\n")
            out_shell.write("{}\n".format(align_cmd))
            out_shell.write("echo \"Starting quantification\"\n")
            if not self.keep_all:
                out_shell.write("{}\n".format(quant_cmd))
                out_shell.write("echo \"Starting cleanup step\"\n")
            out_shell.write("{}\n".format(clean_cmd))
            out_shell.write("echo \"Processing done!\"\n")


        IOMethods.execute_cmd(csv_to_fasta_cmd)

        if not os.path.exists(fasta_file) or os.stat(fasta_file).st_size == 0:
            IOMethods.execute_cmd(index_cmd)

        if not os.path.exists(sam_file) or os.stat(sam_file).st_size == 0:
            IOMethods.execute_cmd(align_cmd)

        if not os.path.exists(bam_file) or os.stat(bam_file).st_size == 0:
            IOMethods.execute_cmd(sam_to_bam_cmd)

        if not os.path.exists(quant_file) or os.stat(quant_file).st_size == 0:
            IOMethods.execute_cmd(quant_cmd)
    
        if not self.keep_all:
            IOMethods.execute_cmd(clean_cmd)

        logging.info("Processing complete for {}".format(self.working_dir))



def add_pipeline_args(parser):
    parser.add_argument(
        "-1",
        "--fq1",
        dest="fq1",
        help="Specify path to Read 1 (R1) FASTQ file"
    )
    parser.add_argument(
        "-2",
        "--fq2",
        dest="fq2",
        help="Specify path to Read 2 (R2) FASTQ file"
    )
    parser.add_argument(
        "-b",
        "--bam_file",
        dest="bam",
        help="Specify path to input BAM file as alternative to FASTQ input"
    )
    parser.add_argument(
        "-s",
        "--sequence_tab",
        dest="seq_tab",
        help="Specify the reference sequences as table with colums name, sequence, and position",
        required=True
    )
    parser.add_argument(
        "-o",
        "--output_folder",
        dest="output_folder",
        help="Specify the folder to save the results into.",
        required=True
    )
    parser.add_argument(
        "-d",
        "--bp_distance",
        dest="bp_distance",
        type=int,
        help="Threshold in base pairs for the required overlap size of reads on both sides of the breakpoint for junction/spanning read counting",
        default=10
    )
    parser.add_argument(
        "--allow_mismatches",
        dest="allow_mismatches",
        action="store_true",
        help="Allow mismatches within the region around the breakpoint determined by the bp_distance parameter"
    )
    parser.add_argument(
        "--interval_mode",
        dest="interval_mode",
        action="store_true",
        help="Specify if interval mode shall be used"
    )
    parser.add_argument(
        "-m",
        "--method",
        dest="method",
        choices=["star", "bowtie2", "bwa"],
        help="Specify alignment software to generate the index",
        default="star"
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="num_threads",
        type=int,
        help="Specify number of threads to use for the alignment",
        default=1
    )
    parser.add_argument(
        "--star_cmd_params",
        dest="star_cmd_params",
        help="Specify STAR commandline parameters to use for the alignment",
        default=""
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--keep_bam",
        dest="keep_alignment",
        help="Do not delete alignment in BAM format during clean up step",
        action="store_true"
    )
    group.add_argument(
        "--keep_all",
        dest="keep_all",
        help="Do not perform clean up step after re-quantification",
        action="store_true"
    )
    parser.set_defaults(func=pipeline_command)


def pipeline_command(args):
    pipe = Pipeline(args.fq1, args.fq2, args.bam, args.seq_tab, args.bp_distance, args.output_folder, args.allow_mismatches, args.interval_mode, args.keep_alignment, args.keep_all)
    pipe.run(args.method, args.num_threads, args.star_cmd_params)
