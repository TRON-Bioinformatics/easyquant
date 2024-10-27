"""
Pipeline module for the overall process.
CSV2FASTA -> Indexing -> Alignment -> Counting
"""

import logging
import os

import bp_quant.io.io_methods as IOMethods
from bp_quant.general.version import VERSION



class Pipeline:
    """
    Creates pipeline object that can be run.
    """

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

        with open(os.path.join(self.working_dir, "run_command.sh"), "w", encoding="utf8") as outf:
            outf.write("#!/bin/sh\n\n")
            outf.write(f"version=\"{VERSION}\"\n")
            outf.write(f"python {os.path.realpath(__file__)} \\\n")
            if self.fq1 and self.fq2:
                outf.write(f"-1 {self.fq1} \\\n")
                outf.write(f"-2 {self.fq2} \\\n")
            elif self.bam:
                outf.write(f"-b {self.bam} \\\n")
            outf.write(f"-s {self.seq_tab} \\\n")
            outf.write(f"-o {self.working_dir} \\\n")
            outf.write(f"-d {self.bp_distance} \\\n")
            outf.write(f"-m {method} \\\n")
            outf.write(f"-t {num_threads} \\\n")
            if self.interval_mode:
                outf.write("--interval_mode \\\n")
            if self.allow_mismatches:
                outf.write("--allow_mismatches \\\n")
            if self.skip_singleton:
                outf.write("--skip_singleton \\\n")
            if align_cmd_params:
                outf.write(align_cmd_params)


        logging.info("Executing bpquant %s", VERSION)
        if self.fq1 and self.fq2:
            logging.info("FQ1=%s", self.fq1)
            logging.info("FQ2=%s", self.fq2)
        elif self.bam:
            logging.info("BAM=%s", self.bam)

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
            clean_up_files.append(f"{cram_file}.crai")
            clean_up_files.append(fasta_file)

        #create folders
        IOMethods.create_folder(align_path)


        if not os.path.exists(num_reads_file) or os.stat(num_reads_file).st_size == 0:
            with open(num_reads_file, "w", encoding="utf8") as outf:
                if self.fq1:
                    outf.write(str(int(IOMethods.get_read_count(self.fq1, "fq"))))
                elif self.bam:
                    outf.write(str(int(IOMethods.get_read_count(self.bam, "bam"))))

        csv_to_fasta_cmd = f"bp_quant csv2fasta \
        --input_csv {self.seq_tab} \
        --output_fasta {fasta_file}"


        index_cmd = f"bp_quant index \
        --input_fasta {fasta_file} \
        --index_dir {genome_path} \
        -t {num_threads} \
        --method {method}"

        align_cmd = ""
        custom_params = ""
        if align_cmd_params:
            custom_params = f"--params '{align_cmd_params}'"
        if self.fq1 and self.fq2:
            align_cmd = f"bp_quant align \
            --fq1 {self.fq1} \
            --fq2 {self.fq2} \
            --index_dir {genome_path} \
            --output_dir {align_path} \
            -t {num_threads} \
            -m {method} \
            {custom_params}"
        elif self.bam:
            align_cmd = f"bp_quant align \
            --bam {self.bam} \
            --index_dir {genome_path} \
            --output_dir {align_path} \
            -t {num_threads} \
            -m {method} \
            {custom_params}"

        if method == "star":
            clean_up_files.extend([
                f"{align_path}/Log.progress.out",
                f"{align_path}/SJ.out.tab",
                f"{align_path}/_STARtmp",
            ])


        sam_to_cram_cmd = f"samtools sort \
            -@ {num_threads} \
            -m 2G \
            --reference {fasta_file} \
            -O CRAM \
            -o {cram_file} {sam_file} \
            && samtools index {cram_file}"

        allow_mismatches_str = ""
        interval_mode_str = ""
        skip_singleton_str = ""
        if self.allow_mismatches:
            allow_mismatches_str = " --allow_mismatches"
        if self.interval_mode:
            interval_mode_str = " --interval_mode"
        if self.skip_singleton:
            skip_singleton_str = " --skip_singleton"

        quant_cmd = f"bp_quant count \
            -i {sam_file} \
            -t {self.seq_tab} \
            -d {self.bp_distance} \
            -o {self.working_dir}\
            {allow_mismatches_str}\
            {interval_mode_str}\
            {skip_singleton_str}"

        clean_up_files_str = " ".join(clean_up_files)
        clean_cmd = f"rm -rf {clean_up_files_str}"

        # define bash script in working directory
        shell_script = os.path.join(self.working_dir, "requant.sh")
        # start to write shell script to execute mapping cmd
        with open(shell_script, "w", encoding="utf8") as out_shell:
            out_shell.write("#!/bin/sh\n\n")
            if self.fq1 and self.fq2:
                out_shell.write(f"fq1_in={self.fq1}\n")
                out_shell.write(f"fq2_in={self.fq2}\n")
            elif self.bam:
                out_shell.write(f"bam_in={self.bam}\n")
            out_shell.write(f"working_dir={self.working_dir}\n")
            out_shell.write("echo \"Starting pipeline...\"\n")
            out_shell.write("echo \"Generating index\"\n")
            out_shell.write(f"{index_cmd}\n")
            out_shell.write("echo \"Starting alignment\"\n")
            out_shell.write(f"{align_cmd}\n")
            out_shell.write("echo \"Starting quantification\"\n")
            out_shell.write(f"{quant_cmd}\n")
            if self.keep_aln or self.keep_all:
                out_shell.write("echo \"Starting SAM->CRAM conversion\"\n")
                out_shell.write(f"{sam_to_cram_cmd}\n")
            else:
                out_shell.write("echo \"Starting cleanup step\"\n")
                out_shell.write(f"{clean_cmd}\n")
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

        logging.info("Processing complete for %s", self.working_dir)



def add_pipeline_args(parser):
    """Parses command line arguments."""
    parser.add_argument(
        "-1",
        "--fq1",
        dest="fq1",
        help="Path to Read 1 (R1) FASTQ file"
    )
    parser.add_argument(
        "-2",
        "--fq2",
        dest="fq2",
        help="Path to Read 2 (R2) FASTQ file"
    )
    parser.add_argument(
        "-b",
        "--bam_file",
        dest="bam",
        help="Path to input BAM file as alternative to FASTQ input"
    )
    parser.add_argument(
        "-s",
        "--sequence_tab",
        dest="seq_tab",
        help="Input sequences as table with columns: `name`, `sequence`, and `position` (tab or `;` separated).",
        required=True
    )
    parser.add_argument(
        "-o",
        "--output_folder",
        dest="output_folder",
        help="Output folder to save the results into.",
        required=True
    )
    parser.add_argument(
        "-d",
        "--bp_distance",
        dest="bp_distance",
        type=int,
        help="Threshold in base pairs for the required overlap size \
            of reads on both sides of the breakpoint for junction/spanning read counting",
        default=10
    )
    parser.add_argument(
        "--allow_mismatches",
        dest="allow_mismatches",
        action="store_true",
        help="Allow mismatches within the region around the breakpoint determined \
            by the bp_distance parameter"
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
        help="Alignment software to generate the index",
        default="star"
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="num_threads",
        type=int,
        help="Number of threads to use for the alignment",
        default=1
    )
    parser.add_argument(
        "--alignment_params",
        dest="align_params",
        help="Custom commandline parameters to use for the alignment",
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
    """Calls pipeline with the parsed command line arguments."""
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
