
import pysam

# parameter --------------------------------------------------------------------
ref = "CLDN18_1"
seq_table_file = "example_data/CLDN18_Context_seq.csv"
# bam_path = "example_out_0.2/star/Aligned.sortedByCoord.out.bam"
bam_path = "bwa_aln/example_rna-seq.bam"
bp_dist = 10
max_mismatch_rate = 0.05



def get_seq_to_pos(seq_table_file):
    """
    Parses the sequence table and returns a dict from name to pos
    """

    # seq_to_pos is a dict from sequence name to the position of interest (breakpoint or junction)

    # initialize dict
    seq_to_pos = {}

    # iterate over all lines and fill dict
    with open(seq_table_file) as in_handle:
        for i, line in enumerate(in_handle):
            # skip header line
            if i >= 1:
                sp = line.strip().split(";")
                name = sp[0]
                pos = int(sp[2])
                seq_to_pos[name] = pos
    return seq_to_pos



def get_read_cash(bam_path):
    """
    Get cash of reads per reference sequence name and per read name (read group)

    Returns a dict of dict of pysam AligndSegment objects
    """

    # oben BAM file
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Initialize empty dict:
    read_cash = {}

    # iterate over all alignments in the BAM file
    for read in bam.fetch():

        # ignore umapped reads
        if read.is_unmapped:
            continue

        # get reference sequence name and read name (read group)
        ref = read.reference_name
        read_name = read.query_name

        # if not existing already, initiallize dict of dict of list
        if not ref in read_cash:
            read_cash[ref] = {}
        if not read_name in read_cash[ref]:
            read_cash[ref][read_name] = []

        # add alignment to dict
        read_cash[ref][read_name].append(read)
    # return dict
    return read_cash



def count_reads(seq_to_pos, cash):
    """
    Count reads and returns a dict seq to count

    returns a dict of reference sequences with a list of read counts:
        0 junc junction reds overlapping the position with at least bp_dist base pairs
        1 span spaning pairs are read pairs of which one mate maps to the left and one to the right of pos
        2 anch maximal overlaping size of junction reads
        3 a number of reads mapping to the region in the left of pos
        4 b number of reads mapping to the region in the right of pos
    """

    # define matche bases for get_aligned_pairs()
    MATCH_BASES = ['A', 'C', 'G', 'T']

    # counts is a dict from sequence names to an array [junc, span]
    # initialize counts
    counts = {}
    for seq_name in seq_to_pos:
        counts[seq_name] = [0, 0, 0, 0, 0]

    # iterate over all reference sequences
    for seq_name in seq_to_pos:

        pos = seq_to_pos[seq_name]

        # initialize anchor
        anchor: int = 0

        # iterate over all reads (read group):
        for read_name in cash[seq_name]:

            # initialize read_name specific indicators
            junc = False
            left = False
            right = False

            # alignments = cash[seq_name][read_name]
            # print(seq_name, read_name, len(alignments))

            # iterate over all alignments:
            for aln in cash[seq_name][read_name]:

                # read overlaps the position of interest
                if aln.reference_start <= pos - bp_dist and aln.reference_end >= pos + bp_dist:

                    print("DEBUG", seq_name, read_name, aln)

                    # test that read maps exact (or no del/ins) in pos +/- bp_dist
                    aln_pairs = aln.get_aligned_pairs(with_seq=True)

                    reg_start = pos - bp_dist
                    reg_end = pos + bp_dist - 1

                    q_start = [q for (q, r, s) in aln_pairs if r == reg_start][0]
                    q_end = [q for (q, r, s) in aln_pairs if r == reg_end][0]

                    # check if boundary postions are aligned and if the stretch length matches on read and reference
                    no_ins_or_del = q_start is not None and \
                                    q_end is not None and \
                                    (q_end - q_start) == (reg_end - reg_start)

                    # test if reference sequence are all capital (means match) acording
                    # to https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_aligned_pairs
                    aln_seq = [s for (q, r, s) in aln_pairs if reg_start <= r and r < reg_end]
                    no_snp = all(s in MATCH_BASES for s in aln_seq)
                    if no_ins_or_del and no_snp:

                        #if seq_name == "CLDN18_2":
                        #    print("DEBUG: CLDN18_2: ", read_name)
                        junc = True
                        anchor = max(anchor, min(aln.reference_end - pos, pos - aln.reference_start))

                # read maps left of the position
                if aln.reference_end <= pos:
                    left = True
                # read maps right of the fusion breakpoint
                if aln.reference_start >= pos:
                    right = True


            # update counters
            if junc: counts[seq_name][0] += 1
            # check if read pairs map to the left and the right part
            if left and right: counts[seq_name][1] += 1

            # update counts of left and right regions
            if left: counts[seq_name][3] += 1
            if right: counts[seq_name][4] += 1

        # add anchor
        counts[seq_name][2] = anchor

    return counts

def write_counts(in_file, counts, out_file):
    """
    Write read counts to output file
    """
    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for i, line in enumerate(in_handle):
                # skip header line
                sp = line.strip().split(";")
                if i == 0:
                    sp_out = ["name", "pos", "junc", "span", "anch", "a", "b"]
                else:
                    # get sequence name from first colum
                    seq_name = sp[0]

                    # get sequence counts
                    seq_counts = counts[seq_name]

                    # drob sequence and construct output columns
                    sp_out = [sp[0], sp[2]] + [str(c) for c in seq_counts]

                out_line = "\t".join(sp_out) + "\n"
                out_handle.write(out_line)


# testing ##############################################################################################################
seq_to_pos = get_seq_to_pos(seq_table_file)
cash = get_read_cash(bam_path)
counts = count_reads(seq_to_pos, cash)
out_file = "CLDN18_Context_seq.csv.counts_new3_star.tsv"
write_counts(seq_table_file, counts, out_file)

# debug
seq_name = "CLDN18_2"
read_name = "J00128:23:H323KBBXX:5:2127:31480:41997"
aln = cash[seq_name][read_name][2]
# aln = cash[seq_name][read_name][1]
