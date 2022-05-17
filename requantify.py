#!/usr/bin/env python

from argparse import ArgumentParser
import csv
import sys

import pysam


def get_seq_to_pos(seq_table_file):
    """
    Parses the sequence table and returns a dict from name to pos
    """

    # seq_to_pos is a dict from sequence name to the position of interest (breakpoint or junction)

    # initialize dict
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


def get_reads_from_bam(bam_path):
    """
    Get cache of reads per reference sequence name and per read name (read group)

    Returns a dict of dict of pysam AligndSegment objects
    """

    # oben BAM file
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Initialize empty dict:
    read_cache = {}

    # iterate over all alignments in the BAM file
    for read in bam.fetch():

        # ignore umapped reads
        #if read.is_unmapped or read.flag > 255:
        if read.is_unmapped:
            continue

        # get reference sequence name and read name (read group)
        ref = read.reference_name
        read_name = read.query_name

        # if not existing already, initiallize dict of dict of list
        if not ref in read_cache:
            read_cache[ref] = {}
        if not read_name in read_cache[ref]:
            read_cache[ref][read_name] = []

        # add alignment to dict
        read_cache[ref][read_name].append((read.reference_start, read.reference_end, read.get_aligned_pairs(with_seq=True)))

    # return dict
    return read_cache

def init_counts(seq_to_pos, interval_mode):
    counts = {}
    cov_perc_dict = {}
    cov_mean_dict = {}
    for seq_name in seq_to_pos:
        if interval_mode:
            counts[seq_name] = {}
            cov_perc_dict[seq_name] = {}
            cov_mean_dict[seq_name] = {}
            for interval_name, ref_start, ref_stop in seq_to_pos[seq_name]:
                # junc, span_read, within, coverage_%, coverage_mean
                counts[seq_name][interval_name] = [0, 0, 0, 0, 0]

                cov_perc_dict[seq_name][interval_name] = {}
                cov_mean_dict[seq_name][interval_name] = {}
                for i in range(ref_start, ref_stop, 1):
                    cov_perc_dict[seq_name][interval_name][i] = 0
                    cov_mean_dict[seq_name][interval_name][i] = 0
        else:
            counts[seq_name] = [0, 0, 0, 0, 0]

    return counts, cov_perc_dict, cov_mean_dict


def count_reads(seq_to_pos, cache, bp_dist, reads_file, interval_mode):
    """
    Count reads and returns a dict seq to count
    returns a dict of reference sequences with a list of read counts:
        0 junc junction reads overlapping the position with at least bp_dist base pairs
        1 span spaning pairs are read pairs of which one mate maps to the left and one to the right of pos
        2 anch maximal overlaping size of junction reads
        3 a number of reads mapping to the region in the left of pos
        4 b number of reads mapping to the region in the right of pos

    if interval_mode is on the output will be as follows:
        0 junc junction reads overlapping the end of the interval with at least bp_dist base pairs
        1 span spaning pairs are read pairs of which each mate maps to a different interval
        2 reads that map completely within the interval
    """

    # Open file for read info
    reads_out = open(reads_file, "w")
    
    # counts is a dict from sequence names to an array [junc, span]
    # initialize counts
    counts, cov_perc_dict, cov_mean_dict = init_counts(seq_to_pos, interval_mode)

    # iterate over all reference sequences
    for seq_name in seq_to_pos:

        if len(seq_to_pos[seq_name]) > 2 and not interval_mode:
            print("Specified too many positions of interest in your input file without using the interval mode!")
            print("Please check your input file or use interval mode!")
            sys.exit(1)

        # test if there are no reads for this sequence
        if seq_name not in cache:
            continue

        left_interval = None
        right_interval = None
        if not interval_mode:
            left_interval = seq_to_pos[seq_name][0][0]
            right_interval = seq_to_pos[seq_name][1][0]

        # initialize anchor
        anchor: int = 0
        # iterate over all reads (read group):
        for read_name in cache[seq_name]:
            if len(cache[seq_name][read_name]) != 2:
                continue
                
            r1_start, r1_stop, r1_pairs = cache[seq_name][read_name][0]
            r2_start, r2_stop, r2_pairs = cache[seq_name][read_name][1]

            # Get read information [junc, within, interval]
            r1_info = classify_read(r1_start, r1_stop, r1_pairs, seq_to_pos[seq_name], bp_dist)
            r2_info = classify_read(r2_start, r2_stop, r2_pairs, seq_to_pos[seq_name], bp_dist)

            anchor = max([r1_info["anchor"], r2_info["anchor"], anchor])
            
            r1_type = ""
            r2_type = ""
            left = False
            right = False
            if r1_info["interval"]:
                interval_name = r1_info["interval"]
                if interval_mode:
                    ref_start, ref_stop = interval_name.split("_")
                    for i in range(max(r1_start, int(ref_start)), min(r1_stop, int(ref_stop)), 1):
                        cov_perc_dict[seq_name][interval_name][i] = 1
                        cov_mean_dict[seq_name][interval_name][i] += 1

                # Check if reads are junction reads
                if r1_info["junc"]:
                    if interval_mode:
                        counts[seq_name][interval_name][0] += 1
                    else:
                        counts[seq_name][0] += 1
                    r1_type = "junc"

                if r1_info["within"]:
                    if interval_mode:
                        counts[seq_name][interval_name][2] += 1
                    else:
                        #print(interval_name, right_interval)
                        if interval_name == left_interval:
                            left = True
                            #counts[seq_name][3] += 1
                        elif interval_name == right_interval:
                            right = True
                            #counts[seq_name][4] += 1
                    r1_type = "within"


            if r2_info["interval"]:
                interval_name = r2_info["interval"]
                if interval_mode:
                    ref_start, ref_stop = interval_name.split("_")
                    for i in range(max(r2_start, int(ref_start)), min(r2_stop, int(ref_stop)), 1):
                        cov_perc_dict[seq_name][interval_name][i] = 1
                        cov_mean_dict[seq_name][interval_name][i] += 1

                if r2_info["junc"]:
                    if interval_mode:
                        counts[seq_name][interval_name][0] += 1
                    else:
                        counts[seq_name][0] += 1                
                    r2_type = "junc"

                if r2_info["within"]:
                    if interval_mode:
                        counts[seq_name][interval_name][2] += 1
                    else:
                        if interval_name == left_interval:
                            left = True
                            #counts[seq_name][3] += 1
                        elif interval_name == right_interval:
                            right = True
                            #counts[seq_name][4] += 1
                    r2_type = "within"


            # Check if r1 and r2 form a spanning pair
            if r1_info["within"] and r2_info["within"] and r1_info["interval"] != r2_info["interval"]:
                if interval_mode:
                    counts[seq_name][r1_info["interval"]][1] += 1
                    counts[seq_name][r2_info["interval"]][1] += 1
                else:
                    counts[seq_name][1] += 1
                r1_type = "span"
                r2_type = "span"

            if left:
                counts[seq_name][3] += 1

            if right:
                counts[seq_name][4] += 1

            if not r1_type:
                r1_type = "softjunc"
            if not r2_type:
                r2_type = "softjunc"
            reads_out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(read_name, "R1", seq_name, r1_start, r1_stop, r1_type))
            reads_out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(read_name, "R2", seq_name, r2_start, r2_stop, r2_type))

        # add anchor
        if not interval_mode:
            counts[seq_name][2] = anchor

    if interval_mode:
        for seq_name in seq_to_pos:
            for interval_name, ref_start, ref_stop in seq_to_pos[seq_name]:
                cov_perc = float(sum(cov_perc_dict[seq_name][interval_name].values()))/(ref_stop - ref_start)
                cov_mean = float(sum(cov_mean_dict[seq_name][interval_name].values()))/(ref_stop - ref_start)

                counts[seq_name][interval_name][3] = cov_perc
                counts[seq_name][interval_name][4] = cov_mean


    reads_out.close()
    return counts


def classify_read(aln_start, aln_stop, aln_pairs, intervals, bp_dist):

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
            if no_ins_or_del and no_snp:
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


def write_counts(in_file, counts, out_file, interval_mode):
    """
    Write read counts to output file
    """
    
    # open input file
    with open(in_file) as csvfile:

        # Auto detect dialect of input file
        dialect = csv.Sniffer().sniff(csvfile.readline(), delimiters=";,\t")
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, dialect=dialect)
        
        # open output file
        with open(out_file, "w") as out_handle:
            
            # write header line
            sp_out = None
            if interval_mode:
                sp_out = ["name", "interval", "overlap_stop", 
                          "span_read", "within_interval", 
                          "coverage_perc", "coverage_mean"]
            else:
                sp_out = ["name", "pos", "junc", "span", "anch", "a", "b"]
            out_line = "\t".join(sp_out) + "\n"
            out_handle.write(out_line)
            
            # Iterate over input file rows
            for row in reader:
                name = row["name"]

                # get sequence counts
                seq_counts = counts[name]

                if interval_mode:
                    for interval_name in seq_counts:
                        # drop sequence and construct output columns
                        sp_out = [name, interval_name] + [str(c) for c in seq_counts[interval_name]]

                        # write as otput line to output file
                        out_line = "\t".join(sp_out) + "\n"
                        out_handle.write(out_line)
                else:
                    position = row["position"]
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
    parser.add_argument('-o', '--output', dest='output', help='Output quantification file in tabular format', default='quantification.tsv')
    parser.add_argument('-d', '--bp_distance', dest='bp_distance', type=int, default=10,
                        help='Distance around postion of interest for junction read counts.')
    parser.add_argument('-r', '--reads_file', dest='reads_file', help='Specify file to store reads information to for later plotting', default='reads.tsv')
    parser.add_argument('--interval-mode', dest='interval_mode', action='store_true', help='Specify if interval mode shall be used')
    args = parser.parse_args()

    # parse sequence names and position of interest
    seq_to_pos = get_seq_to_pos(args.seq_table_file)

    # parse mapped reads into cache
    cache = get_reads_from_bam(bam_path=args.input_bam)

    # count reads
    counts = count_reads(seq_to_pos, cache, args.bp_distance, args.reads_file, args.interval_mode)

    # write to output file
    write_counts(args.seq_table_file, counts, args.output, args.interval_mode)



if __name__ == "__main__":
    main()
