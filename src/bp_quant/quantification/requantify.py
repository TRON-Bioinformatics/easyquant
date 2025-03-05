"""
The quantification module counts all reads in the specified intervals 
and reports them into result files.

@author: Patrick Sorn
"""

import gzip
import logging
from os.path import abspath
from os.path import dirname
from os.path import join as join_path
import sys

# pylint: disable=E0401
import pysam # type: ignore


from bp_quant.validation.alignment_info import get_aligner, get_sorting
from bp_quant.validation.alignment_info import is_chimeric_alignment, is_singleton
from bp_quant.validation.alignment_info import is_valid_alignment, is_singleton_from_read_dict
from bp_quant.quantification.counting import init_count_dict, init_cov_dict, calc_coverage
from bp_quant.quantification.counting import count_reads
from bp_quant.io.file_handler import write_line_to_file
from bp_quant.io.file_handler import save_plot_script
from bp_quant.io.file_handler import save_counts
import bp_quant.io.file_headers as headers
from bp_quant.quantification.read_classification import classify_read, is_spanning_pair
from bp_quant.quantification.read_processing import process_secondary_alignments
from bp_quant.io.seq_table import get_seq_to_pos_dict
from bp_quant.plotting import PLOT_SCRIPT


logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

logger = logging.getLogger(__name__)


class Quantification:
    """Quantification class"""
    def __init__(self, seq_table_file, bam_file, output_path,
                 bp_dist, allow_mismatches, interval_mode, skip_singleton):
        self.bam_file = bam_file
        aln_obj = pysam.AlignmentFile(bam_file, "rb")
        self.aligner = get_aligner(aln_obj)
        self.seq_tab_file = seq_table_file
        self.output_path = abspath(output_path)
        self.quant_file = join_path(output_path, "quantification.tsv")
        self.read_info_file = join_path(output_path, "read_info.tsv.gz")
        self.read_info_out = gzip.open(self.read_info_file, "wb")
        self.read_info_out.write(("\t".join(headers.READ_INFO_HEADER) + "\n").encode())
        self.bp_dist = bp_dist
        self.allow_mismatches = allow_mismatches
        self.interval_mode = interval_mode
        self.skip_singleton = skip_singleton
        self.seq_to_pos = get_seq_to_pos_dict(seq_table_file)
        self.counts = init_count_dict(self.seq_to_pos, self.interval_mode)
        self.cov_dict = init_cov_dict(self.seq_to_pos)
        self.parse_alignment(aln_obj)
        self.read_info_out.close()
        self.write_results()

        sorting = get_sorting(aln_obj)
        if not sorting in ['unsorted', 'queryname']:
            logger.error('Requantification requires queryname sorted SAM file.')
            sys.exit(1)


    def increment_counts(self, seq_name: str, read_counts: dict) -> None:
        """Increments counts for a single read pair."""
        if self.interval_mode:
            for interval_name in read_counts:
                for (key, val) in read_counts[interval_name].items():
                    self.counts[seq_name][interval_name][key] += val
        else:
            for (key, val) in read_counts.items():
                self.counts[seq_name][key] += val


    def update_coverage(self, aln_pairs: list, seq_name: str, intervals: list) -> None:
        """Updates the coverage dictionary with the read info.

        Args:
            aln_pairs (list): Per-base mapping of aligned segment
            seq_name (str): Name of the reference sequence
            interval_name (str): Name of the interval
        """

        for (interval_name, _, _) in intervals:
            ref_start, ref_end = interval_name.split("_")
            for (q, r, _) in aln_pairs:
                if q and r and int(ref_start) <= r < int(ref_end):
                    self.cov_dict[seq_name][interval_name][r] += 1


    def process_reads(self, all_alignments_of_query_name) -> None:
        """
        Process and quantify all alignments of read pair
        """
        secondary_dict = {}
        r1 = None
        r2 = None
        for read in all_alignments_of_query_name:
            # Skip singleton alignments of STAR and bowtie2 if specified by user
            if is_singleton(read) and self.skip_singleton:
                continue
            # Remove chimeric alignments generated by bowtie2 and keep secondary alignments
            quasi_chimeric = self.aligner=="bowtie2" and read.is_secondary
            exclude_chimeric = is_chimeric_alignment(read)

            if not quasi_chimeric and exclude_chimeric:
                continue

            # Bowtie reports alignments in the following order,
            # therefore the following if statement prevents
            # mismatching issues for read pairs where
            # there are secondary alignments:
            # R1
            # R1_secondary_1
            # R1_secondary_2
            # R1_secondary_3
            # R2
            # ...
            # Secondary alignments will be excluded from the matching process
            # and unmapped mates will be simulated as they are not reported
            # within the bowtie2 alignment
            # Gather read information
            qname = read.query_name
            unmapped = read.is_unmapped
            start = -1
            end = -1
            pairs = []
            if not unmapped:
                start = read.reference_start
                end = read.reference_end
                pairs = []
                try:
                    pairs = read.get_aligned_pairs(with_seq=True)
                except AttributeError:
                    logger.error("Not possible to get read information. Skipping...")
                    continue
            read_dict = {
                "query_name": qname,
                "first_in_pair": read.is_read1,
                "unmapped": unmapped,
                "reference_name": read.reference_name,
                "flag": read.flag,
                "cigar": read.cigarstring,
                "start": start,
                "end": end,
                "pairs": pairs
            }
            if read.is_secondary and self.aligner == "bowtie2":
                # Create dictionary for secondary alignments from bowtie2
                # to bring them into a format that can be processed
                if read.reference_name not in secondary_dict:
                    secondary_dict[read.reference_name] = {}
                if read.query_name not in secondary_dict[read.reference_name]:
                    secondary_dict[read.reference_name][read.query_name] = {"R1": [], "R2": []}

                if read_dict["first_in_pair"]:
                    secondary_dict[read.reference_name][read.query_name]["R1"].append(read_dict)
                else:
                    secondary_dict[read.reference_name][read.query_name]["R2"].append(read_dict)
            else:
                # Determine if alignment record comes from R1 or R2
                if read.is_read1:
                    r1 = read_dict
                else:
                    r2 = read_dict

            if r1 and r2:
                # Ignore unmapped read pairs and read pairs
                # with unmatching reference or read name
                if (is_valid_alignment(r1, r2) and
                    is_singleton_from_read_dict(r1, r2)):
                    (seq_name, read_counts) = self.quantify(r1, r2)
                    self.increment_counts(seq_name, read_counts)
                r1 = None
                r2 = None
        # Process secondary alignments individually
        read_pairings = process_secondary_alignments(secondary_dict)
        for r1, r2 in read_pairings:
            # return counts for single read pair and add them up accordingly
            (seq_name, read_counts) = self.quantify(r1, r2)
            self.increment_counts(seq_name, read_counts)


    def parse_alignment(self, aln_handle: pysam.AlignmentFile) -> None:
        """
        Parses alignment and iterates over each read while quantifying it.
        """

        logger.info("Reading alignment file (path=%s).", self.bam_file)
        logger.info("Starting quantification.")

        missing_refs = {}
        all_alignments_of_query_name = []
        current_query_name = None
        for read in aln_handle.fetch():
            # Skip alignments with following SAM flags:
            # * read fails platform/vendor quality checks (0x200)
            # * read is PCR or optical duplicate (0x400)
            # * supplementary alignment (0x800)
            if read.flag > 511:
                continue
            # Handle missing reference sequences which occur in SAM/BAM
            # but not in seq_table.csv
            if read.reference_name and read.reference_name not in self.seq_to_pos:
                if read.reference_name not in missing_refs:
                    missing_refs[read.reference_name] = 0
                missing_refs[read.reference_name] += 1
            # Collect all alignments with the same queryname (same read pair)
            # into chunk and process together before moving to the next read chunk
            if current_query_name is None:
                current_query_name = read.query_name

            if read.query_name == current_query_name:
                all_alignments_of_query_name.append(read)
            else:
                self.process_reads(all_alignments_of_query_name)
                current_query_name = read.query_name
                all_alignments_of_query_name = [read, ]
        # Explicitly call the processing for the last read chunk after the loop
        if all_alignments_of_query_name:
            self.process_reads(all_alignments_of_query_name)

        logger.info("Quantification done.")

        if self.interval_mode:
            self.counts = calc_coverage(self.seq_to_pos, self.cov_dict, self.counts)

        for seq, count in missing_refs.items():
            logger.warning("Could not find %s in reference sequence table: %s reads", seq, count)


    def quantify(self, r1: dict, r2: dict) -> tuple:
        """Quantify read pair r1 and r2 and increment respective count dict values.

        Args:
            r1 (dict): Information on Read 1 of a pair
            r2 (dict): Information on Read 2 of a pair

        Returns:
            tuple: Tuple containing the reference name and the count values
                    for the specific read pair
        """

        seq_name = r1["reference_name"] if not r1["unmapped"] else r2["reference_name"]

        if seq_name not in self.seq_to_pos:
            return seq_name, {}

        intervals = self.seq_to_pos[seq_name][0]

        #logger.info("Read: %s", read_name)
        #logger.info("R1: %s, %s", r1_start, r1_end)
        #logger.info(r1_pairs)
        #logger.info("R2: %s, %s", r2_start, r2_end)
        #logger.info(r2_pairs)

        # Get read information [junc, within, interval]
        r1_info = classify_read(
            aln_pairs=r1["pairs"],
            intervals=intervals,
            bp_dist=self.bp_dist
        )
        r2_info = classify_read(
            aln_pairs=r2["pairs"],
            intervals=intervals,
            bp_dist=self.bp_dist
        )

        if is_spanning_pair(
            r1_info["class"],
            r1_info["interval"],
            r2_info["class"],
            r2_info["interval"]
        ):
            r1_info["class"] = "span"
            r2_info["class"] = "span"

        read_counts = count_reads(
            r1_info, r2_info, intervals, self.allow_mismatches, self.interval_mode
        )

        self.update_coverage(r1["pairs"], seq_name, intervals)
        self.update_coverage(r2["pairs"], seq_name, intervals)

        if r1["unmapped"]:
            r1_info["class"] = "unmapped"
        if r2["unmapped"]:
            r2_info["class"] = "unmapped"

        bp = self.seq_to_pos[seq_name][1]
        write_line_to_file(
            self.read_info_out,
            (
                r1["query_name"],
                "R1", 
                r1["flag"],
                seq_name,
                bp,
                r1["start"],
                r1["end"],
                r1["cigar"],
                r1_info["nm"],
                r1_info["nm_in_bp_area"],
                r1_info["class"],
                r1_info["contains_snp_or_indel_in_bp_area"]
            )
        )
        write_line_to_file(
            self.read_info_out,
            (
                r1["query_name"],
                "R2",
                r2["flag"],
                seq_name,
                bp,
                r2["start"],
                r2["end"],
                r2["cigar"],
                r2_info["nm"],
                r2_info["nm_in_bp_area"],
                r2_info["class"],
                r2_info["contains_snp_or_indel_in_bp_area"]
            )
        )

        return seq_name, read_counts


    def write_results(self) -> None:
        """
        Write read counts and read info to separate output files.
        """

        save_plot_script(
            abspath(self.read_info_file),
            abspath(self.seq_tab_file),
            join_path(self.output_path, "quantification.pdf"),
            join_path(self.output_path, "generate_plot.sh")
        )

        logger.info("Writing results to file (path=%s).", self.quant_file)

        save_counts(self.quant_file, self.seq_to_pos, self.counts, self.interval_mode)


def add_requantify_args(parser):
    """Add arguments for requantification to a parser"""
    parser.add_argument(
        "-i",
        "--input_bam",
        dest="input_bam",
        help="Input BAM file",
        required=True
    )
    parser.add_argument(
        "-t",
        "--seq_table",
        dest="seq_table_file",
        help="Path to input sequence table",
        required=True
    )
    parser.add_argument(
        "-o",
        "--output_path",
        dest="output_path",
        help="Output path where results are stored",
        default="test_out"
    )
    parser.add_argument(
        "-d",
        "--bp_distance",
        dest="bp_distance",
        type=int,
        default=10,
        help="Distance around postion of interest for junction read counts."
    )
    parser.add_argument(
        "--allow_mismatches",
        dest="allow_mismatches",
        action="store_true",
        help="Allow mismatches within the region around the breakpoint \
            determined by the bp_distance parameter"
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
    parser.set_defaults(func=requantify_command)


def requantify_command(args):
    """Run requantification from command line"""
    Quantification(
        seq_table_file=args.seq_table_file,
        bam_file=args.input_bam,
        output_path=args.output_path,
        bp_dist=args.bp_distance,
        allow_mismatches=args.allow_mismatches,
        interval_mode=args.interval_mode,
        skip_singleton=args.skip_singleton
    )
