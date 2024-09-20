"""
Counting module.
"""

import logging
from statistics import mean, median
import sys

from bp_quant.coverage import perc_true

def init_count_dict(seq_to_pos: dict, interval_mode: bool) -> dict:
    """Initializes count dict for further use."""
    counts = {}
    for seq_name in seq_to_pos:
        if interval_mode:
            counts[seq_name] = {}
            for interval_name, _, _ in seq_to_pos[seq_name][0]:
                # junc, span_read, within, coverage_%, coverage_mean, coverage_median
                counts[seq_name][interval_name] = [0, 0, 0, 0, 0, 0]
        else:
            if len(seq_to_pos[seq_name]) > 2:
                logging.error("Specified too many positions of interest \
                    in your input file without using the interval mode!")
                logging.error("Please check your input file or use interval mode!")
                sys.exit(1)

            counts[seq_name] = [0, 0, 0, 0, 0]

    return counts


def init_cov_dict(seq_to_pos: dict) -> dict:
    """Initializes coverage dictionary for further use."""
    cov_dict = {}
    for seq_name in seq_to_pos:
        cov_dict[seq_name] = {}
        for interval_name, ref_start, ref_stop in seq_to_pos[seq_name][0]:

            cov_dict[seq_name][interval_name] = {}
            for i in range(ref_start, ref_stop, 1):
                cov_dict[seq_name][interval_name][i] = 0

    return cov_dict


def calc_coverage(seq_to_pos: dict, cov_dict: dict, counts: dict) -> dict:
    """Calculates the coverage using the values in the coverage dictionary."""
    for seq_name in seq_to_pos:
        seq_vals = seq_to_pos[seq_name]
        for (interval_name, _, _) in seq_vals[0]:
            interval_vals = cov_dict[seq_name][interval_name].values()
            cov_perc = perc_true(interval_vals)
            cov_mean = mean(interval_vals)
            cov_median = median(interval_vals)

            counts[seq_name][interval_name][3] = cov_perc
            counts[seq_name][interval_name][4] = cov_mean
            counts[seq_name][interval_name][5] = cov_median

    return counts


def get_spanning_intervals(intervals: list, r1_interval: str, r2_interval: str) -> list:
    """Returns the intervals the spanning pair uses."""
    start = False
    spanning_intervals = []
    for (interval_name, _, _) in intervals:
        if interval_name == r1_interval or interval_name == r2_interval:
            if not start:
                start = True
            else:
                start = False
        if start:
            spanning_intervals.append(interval_name)
    return spanning_intervals
