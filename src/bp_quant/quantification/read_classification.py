"""
Read classification module.
Reads can be classified into the following types:
junc, within
"""

# define match bases for get_aligned_pairs()
MATCH_BASES = ('A', 'C', 'G', 'T')


def calc_anchor(aln_start: int, aln_end: int, ref_end: int) -> int:
    """Calculates the overlap between alignment and reference stop position."""

    if aln_start < 0:
        raise ValueError(f"aln_start can't be lower than zero ({aln_start})")
    if aln_end < 0:
        raise ValueError(f"aln_end can't be lower than zero ({aln_end})")
    if ref_end < 0:
        raise ValueError(f"ref_end can't be lower than zero ({ref_end})")
    if aln_start > aln_end:
        raise ValueError(f"aln_start cannot be greater than aln_end ({aln_start} > {aln_end})")
    if aln_start == aln_end:
        raise ValueError("aln_start and aln_end cannot be equal ({aln_start} == {aln_end})")

    if aln_start > ref_end:
        return -1
    anchor = max(min(aln_end - ref_end, ref_end - aln_start), 0)
    return anchor


def get_query_pos(pos: int, aln_pairs: list) -> int:
    """Gets the query position on the reference."""
    for (q, r, _) in aln_pairs:
        if r == pos:
            return q
    return -1


def get_match_list(region_start: int, region_end: int, aln_pairs: list) -> list:
    """Extracts the sequence from aln_pairs list within a specified region."""
    aln_seq = []
    if aln_pairs:
        for (_, r, s) in aln_pairs:
            if r and region_start <= r < region_end:
                aln_seq.append(s)
    return aln_seq


def count_mismatches_in_region(region_start: int, region_end: int, aln_pairs: list) -> int:
    """Counts the number of mismatches for a given region."""
    # test if reference sequence are all capital (means match) according
    # to https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_aligned_pairs
    aln_seq = get_match_list(region_start, region_end, aln_pairs)
    match_list = [s in MATCH_BASES for s in aln_seq]
    num_mismatches = match_list.count(False)
    return num_mismatches


def check_for_indels(region_start: int, region_end: int, aln_pairs: list) -> bool:
    """Checks if region contains INDELs."""
    # test that read maps exact (or no del/ins) in pos +/- bp_dist

    q_start = get_query_pos(region_start, aln_pairs)
    q_end = get_query_pos(region_end, aln_pairs)

    # check if boundary positions are aligned
    # and if the stretch length matches on read and reference
    ins_or_del = q_start and \
                    q_end and \
                    (q_end - q_start) != (region_end - region_start)
    return ins_or_del


def get_read_type(bp_dist: int, anchor: int) -> str:
    """Returns the read type based on anchor size."""
    # read overlaps junction area with at least 'bp_dist' bases
    if anchor >= bp_dist:
        return "junc"

    # read overlaps junction area with less than 'bp_dist' bases
    if 0 < anchor < bp_dist:
        return "softjunc"

    # read maps completely within the interval
    if anchor == 0:
        return "within"

    return ""


def classify_read(aln_pairs: list, intervals: list, bp_dist: int) -> dict:
    """Classifies read and returns dict with information on mapping position.
    
    Args:
            aln_pairs (dict): Information on Read 1 of a pair
            r2 (dict): Information on Read 2 of a pair

        Returns:
            tuple: Tuple containing the reference name and the count values
                    for the specific read pair
    """


    read_info = {
        "class": "unclassified", 
        "interval": "", 
        "anchor": 0, 
        "nm": 0, 
        "nm_in_bp_area": 0, 
        "contains_snp_or_indel_in_bp_area": False
    }

    if not aln_pairs:
        return read_info

    aln_start = aln_pairs[0][1]
    aln_end = aln_pairs[-1][1] + 1

    if aln_start < 0:
        return read_info
    # Can we extract start and stop position of alignment
    # from aln_pairs to remove redundancy?
    num_mismatches_total = count_mismatches_in_region(aln_start, aln_end, aln_pairs)

    for (interval_name, _, ref_end) in intervals:

        anchor = calc_anchor(aln_start, aln_end, ref_end)

        read_type = get_read_type(bp_dist, anchor)
        if read_type:
            if read_type == "junc":
                # Get junction area
                reg_start = ref_end - bp_dist
                reg_end = ref_end + bp_dist - 1

                contains_indels = check_for_indels(reg_start, reg_end, aln_pairs)
                num_mismatches_in_bp_area = count_mismatches_in_region(
                    reg_start, reg_end, aln_pairs
                )
                snp_or_indel_in_bp_area = contains_indels or num_mismatches_in_bp_area > 0
                read_info["nm_in_bp_area"] = num_mismatches_in_bp_area
                read_info["contains_snp_or_indel_in_bp_area"] = snp_or_indel_in_bp_area

            read_info["class"] = read_type
            read_info["anchor"] = anchor
            read_info["nm"] = num_mismatches_total
            read_info["interval"] = interval_name
            return read_info

    return read_info


def is_spanning_pair(r1_class: str, r1_interval: str, r2_class: str, r2_interval: str) -> bool:
    """Check if r1 and r2 form a spanning pair."""
    if (r1_class == "within" and
        r2_class == "within" and
        r1_interval != r2_interval
    ):
        return True
    return False
