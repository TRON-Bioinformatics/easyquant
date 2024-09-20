# define match bases for get_aligned_pairs()
MATCH_BASES = ('A', 'C', 'G', 'T')


def calc_anchor(aln_start: int, aln_end: int, ref_end: int):
    """Calculates the overlap between alignment and reference stop position."""

    if aln_start < 0:
        raise ValueError("aln_start can't be lower than zero")
    if aln_end < 0:
        raise ValueError("aln_end can't be lower than zero")
    if ref_end < 0:
        raise ValueError("ref_end can't be lower than zero")
    if aln_start > aln_end:
        raise ValueError("aln_start cannot be greater than aln_end")
    if aln_start == aln_end:
        raise ValueError("aln_start and aln_end cannot be equal")

    if aln_start > ref_end:
        return -1
    anchor = max(min(aln_end - ref_end, ref_end - aln_start), 0)
    return anchor


def get_query_pos(pos: int, aln_pairs: list) -> int:
    """Gets the query position on the reference."""
    for (q, r, s) in aln_pairs:
        if r == pos:
            return q
    return -1
    #return [q for (q, r, s) in aln_pairs if r == pos][0]


def get_match_list(region_start: int, region_end: int, aln_pairs: list) -> list:
    aln_seq = []
    for (q, r, s) in aln_pairs:
        if r and region_start <= r < region_end:
            aln_seq.append(s)
    return aln_seq


def count_mismatches_in_region(region_start: int, region_end: int, aln_pairs: list) -> int:
    
    # test if reference sequence are all capital (means match) according
    # to https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_aligned_pairs
    aln_seq = get_match_list(region_start, region_end, aln_pairs)
    match_list = [s in MATCH_BASES for s in aln_seq]
    num_mismatches = match_list.count(False)
    return num_mismatches


def check_for_indels(region_start: int, region_end: int, aln_pairs: list) -> bool:
    # test that read maps exact (or no del/ins) in pos +/- bp_dist

    q_start = get_query_pos(region_start, aln_pairs)
    q_end = get_query_pos(region_end, aln_pairs)

    # check if boundary positions are aligned and if the stretch length matches on read and reference
    ins_or_del = q_start and \
                    q_end and \
                    (q_end - q_start) != (region_end - region_start)
    return ins_or_del


def get_read_type(bp_dist: int, anchor: int):
    # read overlaps junction area with at least 'bp_dist' bases
    if anchor >= bp_dist:
        return "junc"

    # read overlaps junction area with less than 'bp_dist' bases
    elif 0 < anchor < bp_dist:
        return "softjunc"

    # read maps completely within the interval
    elif anchor == 0:
        return "within"
    else:
        return ""


def classify_read(aln_start, aln_end, aln_pairs, intervals, bp_dist):
    """
    Classifies read and returns dict with information on mapping position
    """

    read_info = {
        "class": "unclassified", 
        "interval": "", 
        "anchor": 0, 
        "nm": 0, 
        "nm_in_bp_area": 0, 
        "contains_snp_or_indel_in_bp_area": False
    }

    # TODO: Can we extract start and stop position of alignment from aln_pairs to remove redundancy?
    # Check performance
    num_mismatches_total = count_mismatches_in_region(aln_start, aln_end, aln_pairs)
    #aln_start = get_query_start_or_end(0, aln_pairs)
    #aln_stop = get_query_start_or_end(0, aln_pairs)

    for (interval_name, ref_start, ref_end) in intervals:
        
        anchor = calc_anchor(aln_start, aln_end, ref_end)

        read_type = get_read_type(bp_dist, anchor)
        if read_type:
            if read_type == "junc":
                # Get junction area
                reg_start = ref_end - bp_dist
                reg_end = ref_end + bp_dist - 1

                contains_indels = check_for_indels(reg_start, reg_end, aln_pairs)
                num_mismatches_in_bp_area = count_mismatches_in_region(reg_start, reg_end, aln_pairs)
                read_info["nm_in_bp_area"] = num_mismatches_in_bp_area
                read_info["contains_snp_or_indel_in_bp_area"] = contains_indels or num_mismatches_in_bp_area > 0

            read_info["class"] = read_type
            read_info["anchor"] = anchor
            read_info["nm"] = num_mismatches_total
            read_info["interval"] = interval_name
            return read_info
        
    return None
    