from argparse import ArgumentParser
import csv
import gzip
import logging
import os
import pysam
from statistics import mean, median
import sys

csv.field_size_limit(sys.maxsize)


logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

logger = logging.getLogger(__name__)


def get_aligner(bam_file):
    """
    Uses pysam to detect the aligner used to create the input BAM file
    """
    header_dict = pysam.AlignmentFile(bam_file, "rb").header.to_dict()
    aligner = header_dict["PG"][0]["ID"].lower()
    return aligner


def perc_true(lst):
    n = len(lst)
    num_true = sum(1 for val in lst if val > 0)
    return float(num_true)/n


def get_seq_to_pos(seq_table_file):
    """
    Parses the sequence table and returns a dict where each sequence name is mapped to intervals
    """

    logger.info("Parsing input sequences from file (path={}).".format(seq_table_file))

    # seq_to_pos is a dict from sequence name to the position of interest (breakpoint or junction)
    seq_to_pos = {}

    with open(seq_table_file, "r") as csvfile:
        # Auto detect dialect of input file
        dialect = csv.Sniffer().sniff(csvfile.readline(), delimiters=";,\t")
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, dialect=dialect)

        # Iterate over input file rows
        for row in reader:
            name = row["name"]

            pos_arr = None
            
            if "," in row["position"]:
                pos_arr = row["position"].split(",")
            else:
                pos_arr = [0, row["position"], len(row["sequence"])]
            intervals = []
            for i in range(len(pos_arr)-1):
                interval_name = "{}_{}".format(pos_arr[i], pos_arr[i+1])
                intervals.append((interval_name, int(pos_arr[i]), int(pos_arr[i+1])))
            seq_to_pos[name] = intervals

    return seq_to_pos



def classify_read(aln_start, aln_stop, aln_pairs, intervals, allow_mismatches, bp_dist):
    """
    Classifies read and returns dict with information on mapping position
    """

    # TODO: Can we extract start and stop position of alignment from aln_pairs to remove redundancy?
    # Check performance
    
    # define match bases for get_aligned_pairs()
    MATCH_BASES = ['A', 'C', 'G', 'T']

    read_info = {"junc": False, "within": False, "interval": "", "anchor": 0, "nm_in_bp_area": 0}

    for (interval_name, ref_start, ref_stop) in intervals:
        # Check if read spans ref start
        if aln_start <= ref_stop - bp_dist and aln_stop >= ref_stop + bp_dist:
            num_mismatches = 0
            # test that read maps exact (or no del/ins) in pos +/- bp_dist

            reg_start = ref_stop - bp_dist
            reg_end = ref_stop + bp_dist - 1

            q_start = [q for (q, r, s) in aln_pairs if r == reg_start][0]
            q_end = [q for (q, r, s) in aln_pairs if r == reg_end][0]

            # check if boundary positions are aligned and if the stretch length matches on read and reference
            no_ins_or_del = q_start is not None and \
                            q_end is not None and \
                            (q_end - q_start) == (reg_end - reg_start)

            # test if reference sequence are all capital (means match) according
            # to https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_aligned_pairs
            aln_seq = [s for (q, r, s) in aln_pairs if r is not None and reg_start <= r and r < reg_end]
            # Get number of mismatches for this read
            match_list = [s in MATCH_BASES for s in aln_seq]
            no_snp = all(match_list)
            num_mismatches = match_list.count(False)
            if (no_ins_or_del and no_snp) or allow_mismatches:
                anchor = min(aln_stop - ref_stop, ref_stop - aln_start)
                read_info["junc"] = True
                read_info["anchor"] = anchor
                read_info["interval"] = interval_name
                read_info["nm_in_bp_area"] = num_mismatches

        
        # read maps completely within the interval
        #if aln_start > ref_start - bp_dist and aln_stop < ref_stop + bp_dist:
        if aln_start >= ref_start and aln_stop <= ref_stop:
            read_info["within"] = True
            read_info["interval"] = interval_name

    return read_info




class Quantification(object):
    def __init__(self, seq_table_file, bam_file, output_path, bp_dist, allow_mismatches, interval_mode):
        self.seq_table_file = seq_table_file
        self.bam_file = bam_file
        self.aligner = get_aligner(bam_file)
        self.output_path = os.path.abspath(output_path)
        self.quant_file = os.path.join(output_path, "quantification.tsv")
        self.reads_file = os.path.join(output_path, "read_info.tsv.gz")
        self.reads_out = gzip.open(self.reads_file, "wb")
        self.reads_out.write(("\t".join([
            "name",
            "mate",
            "reference",
            "start",
            "end",
            "num_mismatches_in_bp_area",
            "classification"
        ]) + "\n").encode())
        self.bp_dist = bp_dist
        self.allow_mismatches = allow_mismatches
        self.interval_mode = interval_mode
        self.seq_to_pos = get_seq_to_pos(seq_table_file)
        self.counts, self.cov_dict = self.init_counts(self.seq_to_pos)
        self.parse_alignment()
        self.reads_out.close()
        self.write_results()


    def init_counts(self, seq_to_pos):
        """
        Initialize count dict.
        """

        logger.info("Initializing count dicts.")
        counts = {}
        cov_dict = {}
        for seq_name in seq_to_pos:

            if self.interval_mode:
                counts[seq_name] = {}
                cov_dict[seq_name] = {}
                for interval_name, ref_start, ref_stop in seq_to_pos[seq_name]:
                    # junc, span_read, within, coverage_%, coverage_mean, coverage_median
                    counts[seq_name][interval_name] = [0, 0, 0, 0, 0, 0]

                    cov_dict[seq_name][interval_name] = {}
                    for i in range(ref_start, ref_stop, 1):
                        cov_dict[seq_name][interval_name][i] = 0
            else:

                if len(seq_to_pos[seq_name]) > 2:
                    logger.error("Specified too many positions of interest in your input file without using the interval mode!")
                    logger.error("Please check your input file or use interval mode!")
                    sys.exit(1)

                counts[seq_name] = [0, 0, 0, 0, 0]

        return counts, cov_dict


    def parse_alignment(self):
        """
        Parses alignment and iterates over each read while quantifying it.
        """

        # TODO: Implement method to parse BAM files using mate information

        logger.info("Reading alignment file (path={}).".format(self.bam_file))
        logger.info("Starting quantification.".format(self.bam_file))
        bam = pysam.AlignmentFile(self.bam_file, "rb")

        missing_refs = {}
        r1 = None
        r2 = None
        for read in bam.fetch():
            if read.flag > 511:
                continue
            # Handle missing reference sequences which occur in SAM/BAM
            # but not in seq_table.csv
            if read.reference_name and read.reference_name not in self.seq_to_pos:
                if read.reference_name not in missing_refs:
                    missing_refs[read.reference_name] = 0
                missing_refs[read.reference_name] += 1
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
            if read.is_secondary and self.aligner == "bowtie2":
                r1_sec = read
                r2_sec = pysam.AlignedSegment(bam.header)
                r2_sec.reference_name = r1_sec.reference_name
                r2_sec.query_name = r1_sec.query_name
                r2_sec.is_unmapped = True
                self.quantify(r1_sec, r2_sec)
            else:
                if not r1:
                    r1 = read
                elif r1 and not r2:
                    r2 = read
                

            if r1 and r2:
                # Ignore unmapped read pairs and read pairs 
                # with unmatching reference or read name
                # TODO: Check if unmapped reference is '*'
                # and therefore needs more relaxed clause
                if not (r1.is_unmapped and r2.is_unmapped or
                    r1.query_name != r2.query_name):
                    if (r1.reference_name == r2.reference_name or
                        r1.is_unmapped and not r2.is_unmapped or
                        not r1.is_unmapped and r2.is_unmapped):
                        self.quantify(r1, r2)
                elif r1.query_name != r2.query_name:
                    print("MISMATCHING QUERY NAME!", r1, r2)
                    break
                r1 = None
                r2 = None
        logger.info("Quantification done.")

        if self.interval_mode:
            for seq_name in self.seq_to_pos:
                for interval_name, ref_start, ref_stop in self.seq_to_pos[seq_name]:
                    cov_perc = perc_true(self.cov_dict[seq_name][interval_name].values())
                    cov_mean = mean(self.cov_dict[seq_name][interval_name].values())
                    cov_median = median(self.cov_dict[seq_name][interval_name].values())

                    self.counts[seq_name][interval_name][3] = cov_perc
                    self.counts[seq_name][interval_name][4] = cov_mean
                    self.counts[seq_name][interval_name][5] = cov_median

        for seq, count in missing_refs.items():
            logger.warning("Could not find {} in reference sequence table: {} reads".format(seq, count))


    def quantify(self, r1, r2):

        """
        Quantify read pair r1 and r2 and increment respective count dict values:
          0 junc junction reads overlapping the position with at least bp_dist base pairs
          1 span spaning pairs are read pairs of which one mate maps to the left and one to the right of reference position
          2 anch maximal overlaping size of junction reads
          3 a number of reads mapping to the region in the left of reference position
          4 b number of reads mapping to the region in the right of reference position

        if interval_mode is on the output will be as follows:
          0 overlap_interval_end_reads reads overlapping the end of the interval with at least bp_dist base pairs
          1 span_interval_end_pairs read pairs spanning the interval end point
          2 within_interval reads that map completely within the interval or softjunctions
        """

        read_name = r1.query_name
        seq_name = r1.reference_name
        
        if seq_name not in self.seq_to_pos:
            return

        left_interval = None
        right_interval = None

        if not self.interval_mode:
            left_interval = self.seq_to_pos[seq_name][0][0]
            right_interval = self.seq_to_pos[seq_name][1][0]


        r1_start = -1
        r1_stop = -1
        r1_pairs = None
        if not r1.is_unmapped:
            r1_start = r1.reference_start
            r1_stop = r1.reference_end
            r1_pairs = None
            try:
                r1_pairs = r1.get_aligned_pairs(with_seq=True)
            except:
                logger.error("Not possible to get read information. Skipping...")
                return

        r2_start = -1
        r2_stop = -1
        r2_pairs = None
        if not r2.is_unmapped:
            r2_start = r2.reference_start
            r2_stop = r2.reference_end
            r2_pairs = None
            try:
                r2_pairs = r2.get_aligned_pairs(with_seq=True)
            except:
                logger.error("Not possible to get read information. Skipping...")
                return



        # Get read information [junc, within, interval]
        r1_info = classify_read(
            aln_start=r1_start,
            aln_stop=r1_stop,
            aln_pairs=r1_pairs,
            intervals=self.seq_to_pos[seq_name],
            allow_mismatches=self.allow_mismatches,
            bp_dist=self.bp_dist
        )
        r2_info = classify_read(
            aln_start=r2_start,
            aln_stop=r2_stop,
            aln_pairs=r2_pairs,
            intervals=self.seq_to_pos[seq_name],
            allow_mismatches=self.allow_mismatches,
            bp_dist=self.bp_dist
        )
        
        if not self.interval_mode:
            self.counts[seq_name][2] = max([r1_info["anchor"], r2_info["anchor"], self.counts[seq_name][2]])
            
        r1_type = ""
        r2_type = ""

        if r1_info["interval"]:
            interval_name = r1_info["interval"]
            if self.interval_mode:
                ref_start, ref_stop = interval_name.split("_")
                for q, r, s in r1_pairs:
                    if q != None and r != None and int(ref_start) <= r < int(ref_stop):
                        self.cov_dict[seq_name][interval_name][r] += 1

            # Check if reads are junction reads
            if r1_info["junc"]:
                if self.interval_mode:
                    self.counts[seq_name][interval_name][0] += 1
                else:
                    self.counts[seq_name][0] += 1
                r1_type = "junc"

            if r1_info["within"]:
                if self.interval_mode:
                    self.counts[seq_name][interval_name][2] += 1
                else:
                    if interval_name == left_interval:
                        self.counts[seq_name][3] += 1
                    elif interval_name == right_interval:
                        self.counts[seq_name][4] += 1
                r1_type = "within"


        if r2_info["interval"]:
            interval_name = r2_info["interval"]
            if self.interval_mode:
                ref_start, ref_stop = interval_name.split("_")
                for q, r, s in r2_pairs:
                    if q != None and r != None and int(ref_start) <= r < int(ref_stop):
                        self.cov_dict[seq_name][interval_name][r] += 1


            if r2_info["junc"]:
                if self.interval_mode:
                    self.counts[seq_name][interval_name][0] += 1
                else:
                    self.counts[seq_name][0] += 1                
                r2_type = "junc"

            if r2_info["within"]:
                if self.interval_mode:
                    self.counts[seq_name][interval_name][2] += 1
                else:
                    if interval_name == left_interval:
                        self.counts[seq_name][3] += 1
                    elif interval_name == right_interval:
                        self.counts[seq_name][4] += 1
                r2_type = "within"


        # Check if r1 and r2 form a spanning pair
        if r1_info["within"] and r2_info["within"] and r1_info["interval"] != r2_info["interval"]:
            if self.interval_mode:

                start = False
                spanning_intervals = []
                for interval_name, ref_start, ref_stop in self.seq_to_pos[seq_name]:
                    if interval_name == r1_info["interval"] or interval_name == r2_info["interval"]:
                        if not start:
                            start = True
                        else:
                            start = False
                    if start:
                        spanning_intervals.append(interval_name)
                for interval_name in spanning_intervals:
                    self.counts[seq_name][interval_name][1] += 1                        
            else:
                self.counts[seq_name][1] += 1
            r1_type = "span"
            r2_type = "span"

        if not r1_type:
            r1_type = "softjunc"
        if not r2_type:
            r2_type = "softjunc"
        r1_mismatches = r1_info["nm_in_bp_area"]
        r2_mismatches = r2_info["nm_in_bp_area"]
        self.reads_out.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                read_name, "R1", seq_name, r1_start, r1_stop, r1_mismatches, r1_type
            ).encode()
        )
        self.reads_out.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                read_name, "R2", seq_name, r2_start, r2_stop, r2_mismatches, r2_type
            ).encode()
        )


    def write_results(self):
        """
        Write read counts to output file
        """

        # Create plotting script
        cmd = "{} -i {} -t {} -o {}".format(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_reads.py"),
            os.path.abspath(self.reads_file),
            os.path.abspath(self.seq_table_file),
            os.path.join(self.output_path, "quantification.pdf")
        )
        with open(os.path.join(self.output_path, "generate_plot.sh"), "w") as outf:
            outf.write("#/bin/bash\n\n{}".format(cmd))


        logger.info("Writing results to file (path={}).".format(self.quant_file))
        
        # open output file
        with open(self.quant_file, "w") as out_handle:
            
            # write header line
            sp_out = None
            if self.interval_mode:
                sp_out = ["name", "interval", "overlap_interval_end_reads", 
                          "span_interval_end_pairs", "within_interval", 
                          "coverage_perc", "coverage_mean", "coverage_median"]
            else:
                sp_out = ["name", "pos", "junc", "span", "anch", "a", "b"]
            out_line = "\t".join(sp_out) + "\n"
            out_handle.write(out_line)
            
            # Iterate over sequence dictionary
            for name in self.seq_to_pos:
                # get sequence counts
                seq_counts = self.counts[name]

                if self.interval_mode:
                    for interval_name in seq_counts:
                        # drop sequence and construct output columns
                        sp_out = [name, interval_name] + [str(c) for c in seq_counts[interval_name]]

                        # write as otput line to output file
                        out_line = "\t".join(sp_out) + "\n"
                        out_handle.write(out_line)
                else:
                    position = str(self.seq_to_pos[name][0][2])
                    # drop sequence and construct output columns
                    sp_out = [name, position] + [str(c) for c in seq_counts]
                    
                    # write as otput line to output file
                    out_line = "\t".join(sp_out) + "\n"
                    out_handle.write(out_line)



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
        help="Allow mismatches within the region around the breakpoint determined by the bp_distance parameter"
    )
    parser.add_argument(
        "--interval_mode",
        dest="interval_mode",
        action="store_true",
        help="Specify if interval mode shall be used"
    )
    parser.set_defaults(func=requantify_command)


def requantify_command(args):
    """Run requantification from command line"""
    requant = Quantification(
        seq_table_file=args.seq_table_file,
        bam_file=args.input_bam,
        output_path=args.output_path,
        bp_dist=args.bp_distance,
        allow_mismatches=args.allow_mismatches,
        interval_mode=args.interval_mode
    )
