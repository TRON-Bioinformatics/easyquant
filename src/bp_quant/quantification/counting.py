"""
Counting module.
"""

import logging
from statistics import mean, median
import sys

# pylint: disable=E0401
from bp_quant.quantification.coverage import perc_true


def init_count_dict_interval_mode(seq_to_pos: dict) -> dict:
    """Initializes count dict for interval mode."""

    counts = {}
    for seq_name in seq_to_pos:
        counts[seq_name] = {}
        for interval_name, _, _ in seq_to_pos[seq_name][0]:
            counts[seq_name][interval_name] = {
                "overlap_interval_end_reads": 0,
                "span_interval_end_pairs": 0,
                "within_interval": 0,
                "coverage_perc": 0.0,
                "coverage_mean": 0.0,
                "coverage_median": 0.0
            }
    return counts


def init_count_dict_regular_mode(seq_to_pos: dict) -> dict:
    """Initializes count dict for regular mode."""

    counts = {}
    for seq_name in seq_to_pos:
        if len(seq_to_pos[seq_name]) > 2:
            logging.error("Specified too many positions of interest \
                in your input file without using the interval mode!")
            logging.error("Please check your input file or use interval mode!")
            sys.exit(1)

        counts[seq_name] = {
            "junc": 0,
            "span": 0,
            "anch": 0,
            "a": 0,
            "b": 0
        }

    return counts


def init_count_dict(seq_to_pos: dict, interval_mode: bool) -> dict:
    """Initializes count dict for regular mode."""

    if interval_mode:
        return init_count_dict_interval_mode(seq_to_pos)
    return init_count_dict_regular_mode(seq_to_pos)


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


def get_spanning_intervals(intervals: list, r1_interval: str, r2_interval: str) -> list:
    """Returns the intervals the spanning pair uses."""
    start = False
    spanning_intervals = []
    for (interval_name, _, _) in intervals:
        if interval_name in (r1_interval, r2_interval):
            start = not start
        if start:
            spanning_intervals.append(interval_name)
    return spanning_intervals


def count_reads_interval_mode(r1_info: dict, r2_info: dict,
                              intervals: list, allow_mm: bool) -> dict:
    """Counts the reads based on their info in interval mode.

    Args:
        r1_info (dict): Information on Read 1 of pair
        r2_info (dict): Information on Read 2 of pair
        intervals (dict): Dict with intervals
        allow_mm (bool): Select if mismatches are allowed

    Returns:
        dict: Read counts for one read pair
    """
    read_counts = {}
    for (interval_name, _, _) in intervals:
        read_counts[interval_name] = {
            "overlap_interval_end_reads": 0,
            "span_interval_end_pairs": 0,
            "within_interval": 0
        }

    if r1_info["interval"]:
        interval_name = r1_info["interval"]
        # Check if reads are junction reads
        if r1_info["class"] == "junc":
            if not r1_info["contains_snp_or_indel_in_bp_area"] or allow_mm:
                read_counts[interval_name]["overlap_interval_end_reads"] += 1

        if r1_info["class"] in ("within", "span"):
            read_counts[interval_name]["within_interval"] += 1

    if r2_info["interval"]:
        interval_name = r2_info["interval"]
        if r2_info["class"] == "junc":
            if not r2_info["contains_snp_or_indel_in_bp_area"] or allow_mm:
                read_counts[interval_name]["overlap_interval_end_reads"] += 1

        if r2_info["class"] in ("within", "span"):
            read_counts[interval_name]["within_interval"] += 1


    # Count spanning pairs individually
    if r1_info["class"] == "span" and r2_info["class"] == "span":
        spanning_intervals = get_spanning_intervals(
            intervals,
            r1_info["interval"],
            r2_info["interval"]
        )
        for interval_name in spanning_intervals:
            read_counts[interval_name]["span_interval_end_pairs"] += 1

    return read_counts


def count_reads_regular_mode(r1_info: dict, r2_info: dict,
                             intervals: list, allow_mm: bool) -> dict:
    """Counts reads based on their info in regular mode.

    Args:
        r1_info (dict): Information on Read 1 of pair
        r2_info (dict): Information on Read 2 of pair
        intervals (dict): Dict with intervals
        allow_mm (bool): Select if mismatches are allowed

    Returns:
        dict: Read counts for one read pair
    """
    left_interval = intervals[0][0]
    right_interval = intervals[1][0]

    read_counts = {
        "junc": 0,
        "span": 0,
        "a": 0,
        "b": 0
    }
    if r1_info["interval"]:
        interval_name = r1_info["interval"]

        # Check if reads are junction reads
        if r1_info["class"] == "junc":
            if not r1_info["contains_snp_or_indel_in_bp_area"] or allow_mm:
                read_counts["junc"] += 1

        if r1_info["class"] in ("within", "span"):
            if interval_name == left_interval:
                read_counts["a"] += 1
            elif interval_name == right_interval:
                read_counts["b"] += 1


    if r2_info["interval"]:
        interval_name = r2_info["interval"]

        if r2_info["class"] == "junc":
            if not r2_info["contains_snp_or_indel_in_bp_area"] or allow_mm:
                read_counts["junc"] += 1

        if r2_info["class"] in ("within", "span"):
            if interval_name == left_interval:
                read_counts["a"] += 1
            elif interval_name == right_interval:
                read_counts["b"] += 1


    # Count spanning pairs individually
    if r1_info["class"] == "span" and r2_info["class"] == "span":
        read_counts["span"] += 1

    return read_counts


def count_reads(r1_info: dict, r2_info: dict,
                intervals: list, allow_mm: bool, interval_mode: bool) -> dict:
    """Counts reads based on their info.

    Args:
        r1_info (dict): Information on Read 1 of pair
        r2_info (dict): Information on Read 2 of pair
        intervals (dict): Dict with intervals
        allow_mm (bool): Select if mismatches are allowed
        interval_mode (bool): Select if interval mode is used

    Returns:
        dict: Read counts for one read pair
    """
    if interval_mode:
        return count_reads_interval_mode(r1_info, r2_info, intervals, allow_mm)
    return count_reads_regular_mode(r1_info, r2_info, intervals, allow_mm)


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
