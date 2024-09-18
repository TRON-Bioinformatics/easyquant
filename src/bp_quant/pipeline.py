import logging
import os
import subprocess
import sys

import bp_quant.io_methods as IOMethods
from bp_quant.version import version


def get_read_count(infile, format="fq"):
    if format == "fq":
        return __get_read_count_fq(infile)
    elif format == "bam":
        return __get_read_count_bam(infile)

    
def __get_read_count_fq(fq_file):
    """Parses input FASTQ to get read count"""
    ps = subprocess.Popen(("zcat", fq_file), stdout=subprocess.PIPE)
    result = subprocess.check_output(("wc", "-l"), stdin=ps.stdout)
    return int(result) / 4


def __get_read_count_bam(bam_file):
    """Parses input BAM to get read count"""
    result = subprocess.check_output(["samtools", "view", "-c", bam_file])
    return int(result)


class Pipeline(object):
    
    def __init__(self, fq1, fq2, bam, seq_tab, bp_distance, working_dir, allow_mismatches, interval_mode, skip_singleton, keep_aln, keep_all):

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
        self.keep_aln = keep_aln
        self.keep_all = keep_all
        self.skip_singleton = skip_singleton

        
    def run(self, method, num_threads, align_cmd_params):
        """This function runs the pipeline on paired-end FASTQ files."""

        script_call = "python {} {}".format(os.path.realpath(__file__), " ".join(sys.argv[1:]))
        
        with open(os.path.join(self.working_dir, "run_command.sh"), "w") as outf:
            outf.write("#!/bin/sh\n\n")
            outf.write("version=\"{}\"\n".format(version))
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
            if self.skip_singleton:
                outf.write("--skip_singleton \\\n")
            if align_cmd_params:
                outf.write(align_cmd_params)


        logging.info("Executing bpquant {}".format(version))
        if self.fq1 and self.fq2:
            logging.info("FQ1={}".format(self.fq1))
            logging.info("FQ2={}".format(self.fq2))
        elif self.bam:
            logging.info("BAM={}".format(self.bam))

        genome_path = os.path.join(self.working_dir, "index")
        IOMethods.create_folder(genome_path)
        align_path = os.path.join(self.working_dir, "alignment")

        fasta_file = os.path.join(self.working_dir, "context.fa")
        sam_file = os.path.join(align_path, "Aligned.out.sam")
        cram_file = os.path.join(align_path, "Aligned.out.cram")
        quant_file = os.path.join(self.working_dir, "quantification.tsv")
        num_reads_file = os.path.join(self.working_dir, "num_reads.txt")
        # Define files to be deleted after successful run
        clean_up_files = [genome_path, sam_file]
        if not self.keep_aln:
            clean_up_files.append(cram_file)
            clean_up_files.append("{}.crai".format(cram_file))
            clean_up_files.append(fasta_file)
        
        #create folders
        IOMethods.create_folder(align_path)


        if not os.path.exists(num_reads_file) or os.stat(num_reads_file).st_size == 0:
            with open(num_reads_file, "w") as outf:
                if self.fq1:
                    outf.write(str(int(get_read_count(self.fq1, "fq"))))
                elif self.bam:
                    outf.write(str(int(get_read_count(self.bam, "bam"))))

        csv_to_fasta_cmd = "bp_quant csv2fasta \
        --input_csv {} \
        --output_fasta {}".format(
            self.seq_tab,
            fasta_file
        )

        index_cmd = "bp_quant index \
        --input_fasta {} \
        --index_dir {} \
        -t {} \
        --method {}".format(
            fasta_file,
            genome_path,
            num_threads,
            method
        )
        
        align_cmd = ""
        custom_params = ""
        if align_cmd_params:
            custom_params = "--params '{}'".format(align_cmd_params)
        if self.fq1 and self.fq2:
            align_cmd = "bp_quant align \
            --fq1 {} \
            --fq2 {} \
            --index_dir {} \
            --output_path {} \
            -t {} \
            -m {} \
            {}".format(
                self.fq1,
                self.fq2,
                genome_path,
                align_path,
                num_threads,
                method,
                custom_params
            )
        elif self.bam:
            align_cmd = "bp_quant align \
            --bam {} \
            --index_dir {} \
            --output_path {} \
            -t {} \
            -m {} \
            {}".format(
                self.bam,
                genome_path,
                align_path,
                num_threads,
                method,
                custom_params
            )

        if method == "star":
            clean_up_files.extend([
                "{}/Log.progress.out".format(align_path),
                "{}/SJ.out.tab".format(align_path),
                "{}/_STARtmp".format(align_path),
            ])

            
        sam_to_cram_cmd = "samtools sort -@ {2} -m 2G --reference {3} -O CRAM -o {1} {0} && samtools index {1}".format(
            sam_file,
            cram_file,
            num_threads,
            fasta_file
        )

        allow_mismatches_str = ""
        interval_mode_str = ""
        skip_singleton_str = ""
        if self.allow_mismatches:
            allow_mismatches_str = " --allow_mismatches"
        if self.interval_mode:
            interval_mode_str = " --interval_mode"
        if self.skip_singleton:
            skip_singleton_str = " --skip_singleton"

        quant_cmd = "bp_quant count \
        -i {0} \
        -t {1} \
        -d {2} \
        -o {3}{4}{5}".format(
            sam_file,
            self.seq_tab,
            self.bp_distance,
            self.working_dir,
            allow_mismatches_str,
            interval_mode_str,
            skip_singleton_str
        )

        clean_cmd = "for file in {}; \
            do rm -rf $file; done".format(" ".join(clean_up_files))

        # define bash script in working directory    
        shell_script = os.path.join(self.working_dir, "requant.sh")
        # start to write shell script to execute mapping cmd
        with open(shell_script, "w") as out_shell:
            out_shell.write("#!/bin/sh\n\n")
            if self.fq1 and self.fq2:
                out_shell.write("fq1_in={}\n".format(self.fq1))
                out_shell.write("fq2_in={}\n".format(self.fq2))
            elif self.bam:
                out_shell.write("bam_in={}\n".format(self.bam))
            out_shell.write("working_dir={}\n".format(self.working_dir))
            out_shell.write("echo \"Starting pipeline...\"\n")
            out_shell.write("echo \"Generating index\"\n")
            out_shell.write("{}\n".format(index_cmd))
            out_shell.write("echo \"Starting alignment\"\n")
            out_shell.write("{}\n".format(align_cmd))
            out_shell.write("echo \"Starting quantification\"\n")
            out_shell.write("{}\n".format(quant_cmd))
            if self.keep_aln or self.keep_all:
                out_shell.write("echo \"Starting SAM->CRAM conversion\"\n")
                out_shell.write("{}\n".format(sam_to_cram_cmd))
            else:
                out_shell.write("echo \"Starting cleanup step\"\n")
                out_shell.write("{}\n".format(clean_cmd))
            out_shell.write("echo \"Processing done!\"\n")


        IOMethods.execute_cmd(csv_to_fasta_cmd)
        IOMethods.execute_cmd(index_cmd)

        if not os.path.exists(sam_file) or os.stat(sam_file).st_size == 0:
            IOMethods.execute_cmd(align_cmd)

        if not os.path.exists(quant_file) or os.stat(quant_file).st_size == 0:
            IOMethods.execute_cmd(quant_cmd)

        if self.keep_aln or self.keep_all:
            if not os.path.exists(cram_file) or os.stat(cram_file).st_size == 0:
                IOMethods.execute_cmd(sam_to_cram_cmd)
            
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
        "--skip_singleton",
        dest="skip_singleton",
        action="store_true",
        help="Skip singleton alignments in requantification"
    )
    parser.add_argument(
        "-m",
        "--method",
        dest="method",
        choices=["star", "bowtie2"],
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
        "--alignment_params",
        dest="align_params",
        help="Specify custom commandline parameters to use for the alignment",
        default=""
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--keep_aln",
        dest="keep_aln",
        help="Do not delete alignment files during clean up step",
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
    pipe = Pipeline(
        fq1=args.fq1,
        fq2=args.fq2,
        bam=args.bam,
        seq_tab=args.seq_tab,
        bp_distance=args.bp_distance,
        working_dir=args.output_folder,
        allow_mismatches=args.allow_mismatches,
        interval_mode=args.interval_mode,
        skip_singleton=args.skip_singleton,
        keep_aln=args.keep_aln,
        keep_all=args.keep_all
    )
    pipe.run(args.method, args.num_threads, args.align_params)
