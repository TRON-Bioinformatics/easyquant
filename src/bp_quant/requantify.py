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

def get_sorting(bam_file):
    """
    Use pysam to grep SO from SAM header and check sorting of input file
    """
    header_dict = pysam.AlignmentFile(bam_file, "rb").header.to_dict()
    sorting = header_dict['HD'].get('SO', 'unsorted')
    return sorting

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
            seq_to_pos[name] = (intervals, ",".join([str(x) for x in pos_arr]))

    return seq_to_pos



def classify_read(aln_start, aln_stop, aln_pairs, intervals, bp_dist):
    """
    Classifies read and returns dict with information on mapping position
    """

    # TODO: Can we extract start and stop position of alignment from aln_pairs to remove redundancy?
    # Check performance
    
    # define match bases for get_aligned_pairs()
    MATCH_BASES = ['A', 'C', 'G', 'T']

    read_info = {"class": "unclassified", "interval": "", "anchor": 0, "nm": 0, "nm_in_bp_area": 0, "contains_snp_or_indel": False}
    for (interval_name, ref_start, ref_stop) in intervals:
        # Check if read spans ref start
        aln_seq_junc = None
        no_ins_or_del = None
        if aln_start <= ref_stop - bp_dist and aln_stop >= ref_stop + bp_dist:
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
            aln_seq_junc = [s for (q, r, s) in aln_pairs if r is not None and reg_start <= r < reg_end]
            read_info["class"] = "junc"
            read_info["interval"] = interval_name
        if aln_start <= ref_stop and aln_stop >= ref_stop:
            aln_seq_junc = [s for (q, r, s) in aln_pairs if r is not None and aln_start <= r < aln_stop]
            # What will happen if multiple BPs with small distance to each other are present?
            # e.g. if a read spans multiple BPs, which anchor will be used?
            # Suggestion: Take the largest anchor among BPs for a read
            # or sum them up
            anchor = min(aln_stop - ref_stop, ref_stop - aln_start)
            if 0 < anchor <= bp_dist:
                read_info["class"] = "softjunc"
            read_info["anchor"] = anchor
            read_info["interval"] = interval_name
        if aln_seq_junc:
            # Get number of mismatches for this read
            match_list_junc = [s in MATCH_BASES for s in aln_seq_junc]
            no_snp = all(match_list_junc)
            num_mismatches_junc = match_list_junc.count(False)
            read_info["nm_in_bp_area"] = num_mismatches_junc
                    # Classify junctions reads independant of INDELs or mismatches
            if no_ins_or_del and no_snp:
                read_info["contains_snp_or_indel"] = False
            else:
                read_info["contains_snp_or_indel"] = True
        num_mismatches = 0
        if aln_pairs:
            aln_seq = [s for (q, r, s) in aln_pairs if r is not None]
            match_list = [s in MATCH_BASES for s in aln_seq]
            num_mismatches = match_list.count(False)
        read_info["nm"] = num_mismatches

        # read maps completely within the interval
        if aln_start >= ref_start and aln_stop <= ref_stop:
            read_info["class"] = "within"
            read_info["interval"] = interval_name
    return read_info


def classify_pair(r1_info, r2_info):
    if r1_info["class"] == "within" and r2_info["class"] == "within" and r1_info["interval"] != r2_info["interval"]:
        r1_info["class"] == "span"
        r2_info["class"] == "span"

    return r1_info, r2_info


def process_secondary_alignments(read_dict):
    """
    Process secondary alignments and classify the reads accordingly
    """

    read_pairings = []
    for key in read_dict:
        for query_name in read_dict[key]:
            # Get every read pairing for the query name
            num_read_pairs = max(
                len(read_dict[key][query_name]["R1"]), 
                len(read_dict[key][query_name]["R2"])
            )
            for i in range(num_read_pairs):
                r1_sec = None
                r2_sec = None
                # Get i-th element of R1 reads, otherwise simulate R2
                try:
                    r1_sec = read_dict[key][query_name]["R1"][i]
                except:
                    r1_sec = {
                        "reference_name": key,
                        "query_name": query_name,
                        "unmapped": True,
                        "flag": 325,
                        "start": -1,
                        "stop": -1,
                        "pairs": None,
                        "cigar": None
                    }
                # Get i-th element of R2 reads, otherwise simulate R2
                try:
                    r2_sec = read_dict[key][query_name]["R2"][i] 
                except:
                    r2_sec = {
                        "reference_name": key,
                        "query_name": query_name,
                        "unmapped": True,
                        "flag": 389,
                        "start": -1,
                        "stop": -1,
                        "pairs": None,
                        "cigar": None
                    }
                if r1_sec and r2_sec:
                    read_pairings.append((r1_sec, r2_sec))
    return read_pairings

def is_chimeric_alignment(read: pysam.AlignedSegment) -> bool:
    """
    Determine if a given alignment is chimeric with mates mapping to different context sequences.
    Following cases have to be considered.

    * read is mapped (reference) != rnext is mapped (reference) -> TRUE
    * read is unmapped (*) and rnext is mapped -> FALSE
    * reference is mapped and rnext is unmapped (*) -> FALSE
    """
    read_reference = read.reference_id
    rnext_reference = read.next_reference_id

    read_mapped = read.is_mapped
    rnext_mapped = read.mate_is_mapped

    if read_mapped and rnext_mapped and read_reference != rnext_reference:
        return True

    return False

def is_singleton(read: pysam.AlignedSegment) -> bool:
    """
    Check if read is a singleton alignment
    """
    if read.is_mapped and read.mate_is_unmapped:
        return True
    elif read.is_unmapped and read.mate_is_mapped:
        return True
    else:
        return False


class Quantification(object):
    def __init__(self, seq_table_file, bam_file, output_path, bp_dist, allow_mismatches, interval_mode, skip_singleton):
        self.seq_table_file = seq_table_file
        self.bam_file = bam_file
        self.aligner = get_aligner(bam_file)
        self.sorting = get_sorting(bam_file)
        self.output_path = os.path.abspath(output_path)
        self.quant_file = os.path.join(output_path, "quantification.tsv")
        self.reads_file = os.path.join(output_path, "read_info.tsv.gz")
        self.reads_out = gzip.open(self.reads_file, "wb")
        self.reads_out.write(("\t".join([
            "name",
            "mate",
            "flag",
            "reference",
            "bp",
            "start",
            "end",
            "cigar",
            "num_mismatches",
            "num_mismatches_in_bp_area",
            "classification",
            "contains_snp_or_indel"
        ]) + "\n").encode())
        self.bp_dist = bp_dist
        self.allow_mismatches = allow_mismatches
        self.interval_mode = interval_mode
        self.skip_singleton = skip_singleton
        self.seq_to_pos = get_seq_to_pos(seq_table_file)
        self.counts, self.cov_dict = self.init_counts(self.seq_to_pos)
        self.parse_alignment()
        self.reads_out.close()
        self.write_results()

        if not self.sorting in ['unsorted', 'queryname']:
            logger.error('Requantification requires queryname sorted SAM file.')
            sys.exit(1)

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
                for interval_name, ref_start, ref_stop in seq_to_pos[seq_name][0]:
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

    def process_reads(self, all_alignments_of_query_name):
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
            # Remove chimeric alignments generated by bowtie2
            if is_chimeric_alignment(read):
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
            stop = -1
            pairs = None
            if not unmapped:
                start = read.reference_start
                stop = read.reference_end
                pairs = None
                try:
                    pairs = read.get_aligned_pairs(with_seq=True)
                except:
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
                "stop": stop,
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
                # TODO: Check if unmapped reference is '*'
                # and therefore needs more relaxed clause
                if not (r1["unmapped"] and r2["unmapped"] or
                    r1["query_name"] != r2["query_name"]):
                    if (r1["reference_name"] == r2["reference_name"] or
                        r1["unmapped"] and not r2["unmapped"] or
                        not r1["unmapped"] and r2["unmapped"]):
                        self.quantify(r1, r2)
                elif r1["query_name"] != r2["query_name"]:
                    print("MISMATCHING QUERY NAME!", r1, r2)
                    break
                r1 = None
                r2 = None
        # Process secondary alignments individually
        read_pairings = process_secondary_alignments(secondary_dict)
        for r1, r2 in read_pairings:
            self.quantify(r1, r2)


    def parse_alignment(self):
        """
        Parses alignment and iterates over each read while quantifying it.
        """

        # TODO: Implement method to parse BAM files using mate information

        logger.info("Reading alignment file (path={}).".format(self.bam_file))
        logger.info("Starting quantification.".format(self.bam_file))
        bam = pysam.AlignmentFile(self.bam_file, "rb")

        missing_refs = {}
        all_alignments_of_query_name = []
        current_query_name = None
        for read in bam.fetch():
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
            for seq_name in self.seq_to_pos:
                for interval_name, ref_start, ref_stop in self.seq_to_pos[seq_name][0]:
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

        read_name = r1["query_name"]
        r1_flag = r1["flag"]
        r2_flag = r2["flag"]
        
        seq_name = r1["reference_name"] if not r1["unmapped"] else r2["reference_name"]
        
        if seq_name not in self.seq_to_pos:
            return

        left_interval = None
        right_interval = None

        if not self.interval_mode:
            left_interval = self.seq_to_pos[seq_name][0][0][0]
            right_interval = self.seq_to_pos[seq_name][0][1][0]


        r1_start = r1["start"]
        r1_stop = r1["stop"]
        r1_pairs = r1["pairs"]

        r2_start = r2["start"]
        r2_stop = r2["stop"]
        r2_pairs = r2["pairs"]

        r1_cigar = r1["cigar"]
        r2_cigar = r2["cigar"]

        # Get read information [junc, within, interval]
        # r1_info = classify_read(
        #     aln_start=r1_start,
        #     aln_stop=r1_stop,
        #     aln_pairs=r1_pairs,
        #     intervals=self.seq_to_pos[seq_name][0],
        #     bp_dist=self.bp_dist
        # )
        # r2_info = classify_read(
        #     aln_start=r2_start,
        #     aln_stop=r2_stop,
        #     aln_pairs=r2_pairs,
        #     intervals=self.seq_to_pos[seq_name][0],
        #     bp_dist=self.bp_dist
        # )

        r1_info, r2_info = classify_pair(classify_read(
            aln_start=r1_start,
            aln_stop=r1_stop,
            aln_pairs=r1_pairs,
            intervals=self.seq_to_pos[seq_name][0],
            bp_dist=self.bp_dist
        ), classify_read(
            aln_start=r2_start,
            aln_stop=r2_stop,
            aln_pairs=r2_pairs,
            intervals=self.seq_to_pos[seq_name][0],
            bp_dist=self.bp_dist
        ))
        
        if not self.interval_mode:
            self.counts[seq_name][2] = max([r1_info["anchor"], r2_info["anchor"], self.counts[seq_name][2]])

        if r1_info["interval"]:
            interval_name = r1_info["interval"]
            if self.interval_mode:
                ref_start, ref_stop = interval_name.split("_")
                for q, r, s in r1_pairs:
                    if q != None and r != None and int(ref_start) <= r < int(ref_stop):
                        self.cov_dict[seq_name][interval_name][r] += 1

            # Check if reads are junction reads
            if r1_info["class"] == "junc":
                if not r1_info["contains_snp_or_indel"] or self.allow_mismatches:
                    if self.interval_mode:
                        self.counts[seq_name][interval_name][0] += 1
                    else:
                        self.counts[seq_name][0] += 1

            if r1_info["class"] in ("within", "span"):
                if self.interval_mode:
                    self.counts[seq_name][interval_name][2] += 1
                else:
                    if interval_name == left_interval:
                        self.counts[seq_name][3] += 1
                    elif interval_name == right_interval:
                        self.counts[seq_name][4] += 1


        if r2_info["interval"]:
            interval_name = r2_info["interval"]
            if self.interval_mode:
                ref_start, ref_stop = interval_name.split("_")
                for q, r, s in r2_pairs:
                    if q != None and r != None and int(ref_start) <= r < int(ref_stop):
                        self.cov_dict[seq_name][interval_name][r] += 1


            if r2_info["class"] == "junc":
                if not r2_info["contains_snp_or_indel"] or self.allow_mismatches:
                    if self.interval_mode:
                        self.counts[seq_name][interval_name][0] += 1
                    else:
                        self.counts[seq_name][0] += 1                

            if r2_info["class"] in ("within", "span"):
                if self.interval_mode:
                    self.counts[seq_name][interval_name][2] += 1
                else:
                    if interval_name == left_interval:
                        self.counts[seq_name][3] += 1
                    elif interval_name == right_interval:
                        self.counts[seq_name][4] += 1


        # Count spanning pairs individually
        if r1_info["class"] == "span" and r2_info["class"] == "span":
            if self.interval_mode:
                start = False
                spanning_intervals = []
                for interval_name, ref_start, ref_stop in self.seq_to_pos[seq_name][0]:
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

        if r1["unmapped"]:
            r1_info["class"] = "unmapped"
        if r2["unmapped"]:
            r2_info["class"] = "unmapped"

        bp = self.seq_to_pos[seq_name][1]

        self.reads_out.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                read_name, "R1", r1_flag, 
                seq_name, bp, 
                r1_start, r1_stop, 
                r1_cigar, 
                r1_info["nm"], 
                r1_info["nm_in_bp_area"], 
                r1_info["class"],
                r1_info["contains_snp_or_indel"]
            ).encode()
        )
        self.reads_out.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                read_name, "R2", r2_flag, 
                seq_name, bp, 
                r2_start, r2_stop, 
                r2_cigar, 
                r2_info["nm"], 
                r2_info["nm_in_bp_area"], 
                r2_info["class"],
                r2_info["contains_snp_or_indel"]
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
                    position = str(self.seq_to_pos[name][0][0][2])
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
    parser.add_argument(
        "--skip_singleton",
        dest="skip_singleton",
        action="store_true",
        help="Skip singleton alignments in requantification"
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
        interval_mode=args.interval_mode,
        skip_singleton=args.skip_singleton
    )
