#!/usr/bin/env python

from argparse import ArgumentParser
import csv

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib.path import Path


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


def count_reads(seq_to_pos, cache, bp_dist):
    """
    Count reads and returns a dict seq to count

    returns a dict of reference sequences with a list of read counts:
        0 junc junction reads overlapping the end of the interval with at least bp_dist base pairs
        1 span spaning pairs are read pairs of which each mate maps to a different interval
        2 reads that map completely within the interval
    """


    # counts is a dict from sequence names to an array [junc, span]
    # initialize counts
    counts = {}
    # Initialize coverage dictionaries
    # coverage calculation will be based on a boolean dict 
    # for percentual coverage counting 1s
    cov_perc_dict = {}
    # Coverage calculation will be based on an integer dict 
    # for fold coverage counting reads overlapping each position
    cov_mean_dict = {}

    for seq_name in seq_to_pos:
        counts[seq_name] = {}
        cov_perc_dict[seq_name] = {}
        cov_mean_dict[seq_name] = {}
        for interval_name, ref_start, ref_stop in seq_to_pos[seq_name]:
            # junc_start, junc_end, span_read, within, coverage_%, coverage_mean
            counts[seq_name][interval_name] = [0, 0, 0, 0, 0]

            cov_perc_dict[seq_name][interval_name] = {}
            cov_mean_dict[seq_name][interval_name] = {}
            for i in range(ref_start, ref_stop, 1):
                cov_perc_dict[seq_name][interval_name][i] = 0
                cov_mean_dict[seq_name][interval_name][i] = 0


    
    pdf = PdfPages('quantification.pdf')
    for seq_name in seq_to_pos:
        if seq_name not in cache:
            continue
        
        fig, ax = plt.subplots()
        #fig, ax = plt.figure(figsize=(3, 3))
        
        plt.axvline(x=200, color='grey', linewidth=1)
        plt.axvline(x=800, color='grey', linewidth=1)

        h = 10
        for read_name in cache[seq_name]:
            r1_start, r1_stop, r1_pairs = cache[seq_name][read_name][0]
            r2_start, r2_stop, r2_pairs = cache[seq_name][read_name][1]
            


            # iterate over all alignments:
            #for start, stop, pairs in cache[seq_name][read_name]:
            #print(r1_start, r1_stop)
            #print(r2_start, r2_stop)
            read_dict = {}
            # Check in all intervals for an overlap
            for interval_name, ref_start, ref_stop in seq_to_pos[seq_name]:
                
                r1_class = classify_read(r1_start, r1_stop, r1_pairs, ref_start, ref_stop, bp_dist)
                r2_class = classify_read(r2_start, r2_stop, r2_pairs, ref_start, ref_stop, bp_dist)
                # read overlaps the position of interest
                #if not r1_class and not r2_class:
                #    continue

                if sum(r1_class):
                    read_dict["1"] = (interval_name, r1_class)

                    for i in range(max(r1_start, ref_start), min(r1_stop, ref_stop), 1):
                        cov_perc_dict[seq_name][interval_name][i] = 1
                        cov_mean_dict[seq_name][interval_name][i] += 1

                    if r1_class[0]:
                        counts[seq_name][interval_name][0] += 1

                    if r1_class[1]:
                        counts[seq_name][interval_name][2] += 1

                if sum(r2_class):
                    read_dict["2"] = (interval_name, r2_class)

                    for i in range(max(r2_start, ref_start), min(r2_stop, ref_stop), 1):
                        cov_perc_dict[seq_name][interval_name][i] = 1
                        cov_mean_dict[seq_name][interval_name][i] += 1
                    # check if read overlaps starting point of interval
                
                    if r2_class[0]:
                        counts[seq_name][interval_name][0] += 1
                    # check if read overlaps stopping point of interval

                    if r2_class[1]:
                        counts[seq_name][interval_name][2] += 1

            


                #print("{}: interval={}, start={}, stop={}, read_class={}".format(read_name, interval_name, start, stop, read_class))
                # update counters

                # check if read maps completely within interval boundaries
                #if read_class == 3:
                #    counts[seq_name][interval_name][2] += 1
            #print(read_dict)
            r1_interval = read_dict["1"][0]
            r2_interval = read_dict["2"][0]
            r1_class = read_dict["1"][1]
            r2_class = read_dict["2"][1]
            r1 = None
            r2 = None
            if r1_interval != r2_interval:
                counts[seq_name][r1_interval][1] += 1
                counts[seq_name][r2_interval][1] += 1
                r1 = Rectangle((r1_start, h), r1_stop-r1_start, 4, color='yellow')
                r2 = Rectangle((r2_start, h), r2_stop-r2_start, 4, color='yellow')
            
            if r1_class[0]:
                r1 = Rectangle((r1_start, h), r1_stop-r1_start, 4, color='red')
            elif r1_class[1] and r1_interval == r2_interval:
                r1 = Rectangle((r1_start, h), r1_stop-r1_start, 4, color='green')
            if r2_class[0]:
                r2 = Rectangle((r2_start, h), r2_stop-r2_start, 4, color='red')
            elif r2_class[1] and r1_interval == r2_interval:
                r2 = Rectangle((r2_start, h), r2_stop-r2_start, 4, color='green')

            plt.plot([r1_stop+2, r2_start-3], [h+2, h+2], color='grey', linestyle='dashed', linewidth=0.5)

            ax.add_patch(r1)
            ax.add_patch(r2)

            h += 12


            # add anchor
            #counts[seq_name][2] = anchor

        plt.xticks(range(0, 1100, 100))
        #plt.yticks(range(0, h+10, 10))
        #plt.show()
        plt.title(seq_name)
        #fig = plt.show()
        pdf.savefig(fig)
        plt.close()
    pdf.close()

    for seq_name in seq_to_pos:
        for interval_name, ref_start, ref_stop in seq_to_pos[seq_name]:
            cov_perc = float(sum(cov_perc_dict[seq_name][interval_name].values()))/(ref_stop - ref_start)
            cov_mean = float(sum(cov_mean_dict[seq_name][interval_name].values()))/(ref_stop - ref_start)

            counts[seq_name][interval_name][3] = cov_perc
            counts[seq_name][interval_name][4] = cov_mean
    print(counts)



    return counts


def classify_read(aln_start, aln_stop, aln_pairs, ref_start, ref_stop, bp_dist):

    # define match bases for get_aligned_pairs()
    MATCH_BASES = ['A', 'C', 'G', 'T']

    read_class = []
    # Check if read spans ref start
    if aln_start <= ref_stop - bp_dist and aln_stop >= ref_stop + bp_dist:
        

        # print("DEBUG", seq_name, read_name, aln)
        
        # test that read maps exact (or no del/ins) in pos +/- bp_dist

        reg_start = ref_stop - bp_dist
        reg_end = ref_stop + bp_dist - 1
        #print(reg_start, reg_end)
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
            read_class.append(1)
        else:
            read_class.append(0)
                #anchor = max(anchor, min(aln_stop - pos, pos - aln_start))
    else:
        read_class.append(0)
        
    # read maps completely within the interval
    if aln_start > ref_start - bp_dist and aln_stop < ref_stop + bp_dist:
        read_class.append(1)
    else:
        read_class.append(0)
    # read maps completely outside the interval
    #elif aln_start < ref_start and aln_stop < ref_start or aln_start > ref_stop and aln_stop > ref_stop:
    #    read_class = 4
    return read_class 


def write_counts(in_file, counts, out_file):
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
            sp_out = ["name", "interval", "overlap_stop", "span_read", "within_interval", "coverage_perc", "coverage_mean"]
            out_line = "\t".join(sp_out) + "\n"
            out_handle.write(out_line)

            # Iterate over input file rows
            for row in reader:

                name = row["name"]
                # sequence = row["sequence"]
                #position = row["position"]

                # get sequence counts
                seq_counts = counts[name]

                for interval_name in seq_counts:

                    # drob sequence and construct output columns
                    sp_out = [name, interval_name] + [str(c) for c in seq_counts[interval_name]]

                    # write as otput line to output file
                    out_line = "\t".join(sp_out) + "\n"
                    out_handle.write(out_line)

    


def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input_bam', dest='input_bam', help='Input BAM file', required=True)
    parser.add_argument('-t', '--seq_table', dest='seq_table_file', help='Path to input sequence table', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Output file', required=True)
    parser.add_argument('-d', '--bp_distance', dest='bp_distance', type=int, default=10,
                        help='Distance around postion of interest for junction read counts.')
    args = parser.parse_args()

    # parse sequence names and position of interest
    seq_to_pos = get_seq_to_pos(args.seq_table_file)

    # parse mapped reads into cache
    cache = get_reads_from_bam(bam_path=args.input_bam)

    # count reads
    counts = count_reads(seq_to_pos, cache, bp_dist=args.bp_distance)
    #print(counts)
    # write to output file
    write_counts(args.seq_table_file, counts, args.output)



if __name__ == "__main__":
    main()
