#!/usr/bin/env python

from argparse import ArgumentParser
import csv
import logging
import os
import sys

import pysam

csv.field_size_limit(sys.maxsize)


logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

logger = logging.getLogger(__name__)


def perc_true(lst):
    n = len(lst)
    num_true = sum(1 for val in lst if val > 0)
    return float(num_true)/n


def mean(lst):
    n = len(lst)
    return float(sum(lst))/n


def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (s[n//2-1]/2.0+s[n//2]/2.0, float(s[n//2]))[n % 2] if n else None


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

    # define match bases for get_aligned_pairs()
    MATCH_BASES = ['A', 'C', 'G', 'T']

    read_info = {"junc": False, "within": False, "interval": "", "anchor": 0}

    for interval_name, ref_start, ref_stop in intervals:
        # Check if read spans ref start
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
            aln_seq = [s for (q, r, s) in aln_pairs if r is not None and reg_start <= r and r < reg_end]
            no_snp = all(s in MATCH_BASES for s in aln_seq)
            if (no_ins_or_del and no_snp) or allow_mismatches:
                anchor = min(aln_stop - ref_stop, ref_stop - aln_start)
                read_info["junc"] = True
                read_info["anchor"] = anchor
                read_info["interval"] = interval_name

        
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
        self.output_path = os.path.abspath(output_path)
        self.quant_file = os.path.join(output_path, "quantification.tsv")
        self.reads_file = os.path.join(output_path, "read_info.tsv")
        self.reads_out = open(self.reads_file, "w")
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

        logger.info("Reading alignment file (path={}).".format(self.bam_file))
        logger.info("Starting quantification.".format(self.bam_file))
        bam = pysam.AlignmentFile(self.bam_file, "rb")

        r1 = None
        r2 = None

        for read in bam.fetch():

            if not r1:
                r1 = read
            elif r1 and not r2:
                r2 = read
                

            if r1 and r2:
                # Ignore unmapped reads and read pairs 
                # with unmatching reference or read name
                if not (r1.is_unmapped or r2.is_unmapped or
                    r1.query_name != r2.query_name or
                    r1.reference_name != r2.reference_name):
                    self.quantify(r1, r2)
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

        left_interval = None
        right_interval = None

        if not self.interval_mode:
            left_interval = self.seq_to_pos[seq_name][0][0]
            right_interval = self.seq_to_pos[seq_name][1][0]


        r1_start = r1.reference_start
        r1_stop = r1.reference_end
        r1_pairs = r1.get_aligned_pairs(with_seq=True)
        r2_start = r2.reference_start
        r2_stop = r2.reference_end
        r2_pairs = r2.get_aligned_pairs(with_seq=True)

        # Get read information [junc, within, interval]
        r1_info = classify_read(r1_start, r1_stop, r1_pairs, self.seq_to_pos[seq_name], self.allow_mismatches, self.bp_dist)
        r2_info = classify_read(r2_start, r2_stop, r2_pairs, self.seq_to_pos[seq_name], self.allow_mismatches, self.bp_dist)

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
        self.reads_out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(read_name, "R1", seq_name, r1_start, r1_stop, r1_type))
        self.reads_out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(read_name, "R2", seq_name, r2_start, r2_stop, r2_type))


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

    

def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input_bam', dest='input_bam', help='Input BAM file', required=True)
    parser.add_argument('-t', '--seq_table', dest='seq_table_file', help='Path to input sequence table', required=True)
    parser.add_argument('-o', '--output_path', dest='output_path', help='Output path where results are stored', default="test_out")
    parser.add_argument('-d', '--bp_distance', dest='bp_distance', type=int, default=10,
                        help='Distance around postion of interest for junction read counts.')
    parser.add_argument('--allow_mismatches', dest='allow_mismatches', action='store_true', help='Allow mismatches within the region around the breakpoint determined by the bp_distance parameter')
    parser.add_argument('--interval_mode', dest='interval_mode', action='store_true', help='Specify if interval mode shall be used')
    args = parser.parse_args()


    q = Quantification(args.seq_table_file, args.input_bam, args.output_path, args.bp_distance, args.allow_mismatches, args.interval_mode)
    

if __name__ == "__main__":
    main()
