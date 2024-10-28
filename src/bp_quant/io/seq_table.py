"""
Module for handling the import of a reference table and generating intervals
from it for downstream counting.
"""

import csv
import sys

csv.field_size_limit(sys.maxsize)


def generate_intervals(pos_arr: list) -> list:
    """Generates intervals from a positions list."""
    intervals = []
    for i in range(len(pos_arr)-1):
        interval_name = f"{pos_arr[i]}_{pos_arr[i+1]}"
        intervals.append((interval_name, int(pos_arr[i]), int(pos_arr[i+1])))
    return intervals


def read_table(seq_table_file: str):
    """Parses input table and yields iterator."""

    with open(seq_table_file, "r", encoding="utf8") as csvfile:
        # Auto detect dialect of input file
        dialect = csv.Sniffer().sniff(csvfile.readline(), delimiters=";,\t")
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, dialect=dialect)
        for row in reader:
            yield row


def get_seq_to_pos_dict(seq_table_file: str) -> dict:
    """
    Parses the sequence table and returns a dict where each sequence name is mapped to intervals
    """

    # seq_to_pos is a dict from sequence name to the position of interest (breakpoint or junction)
    seq_to_pos = {}

    # Iterate over input file rows
    for row in read_table(seq_table_file):
        name = row["name"]

        pos_arr = None

        if "," in row["position"]:
            pos_arr = row["position"].split(",")
        else:
            pos_arr = [0, row["position"], len(row["sequence"])]
        intervals = generate_intervals(pos_arr)
        seq_to_pos[name] = (intervals, ",".join([str(x) for x in pos_arr]))
    return seq_to_pos
