"""
Library for file headers.
"""

COUNTS_INT_MODE = [
    "name", 
    "interval", 
    "overlap_interval_end_reads", 
    "span_interval_end_pairs", 
    "within_interval", 
    "coverage_perc", 
    "coverage_mean", 
    "coverage_median"
]

COUNTS_SINGLE_MODE = [
    "name", 
    "pos", 
    "junc", 
    "span", 
    "anch", 
    "a", 
    "b"
]

READ_INFO_HEADER = [
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
    "contains_snp_or_indel_in_bp_area"
]
