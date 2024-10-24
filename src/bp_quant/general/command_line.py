"""
This module generates the command line subparsers for
the other modules.
"""

import argparse

from bp_quant.io.csv_to_fasta import add_csv_to_fasta_args
from bp_quant.indexing.generate_index import add_indexing_args
from bp_quant.alignment.generate_alignment import add_aligner_args
from bp_quant.plotting.plot_reads import add_plot_reads_args
from bp_quant.general.pipeline import add_pipeline_args
from bp_quant.general.version import VERSION
from bp_quant.quantification.requantify import add_requantify_args


EPILOG = "Copyright (c) 2024 TRON gGmbH (See LICENSE for licensing details)"


def bp_quant_cli():
    """
    This method generates all the parsers.
    """
    parser = argparse.ArgumentParser(
        description=f"bpquant version {VERSION}",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=EPILOG,
    )

    subparsers = parser.add_subparsers(description="Commands")

    pipeline_parser = subparsers.add_parser(
        "pipeline",
        description="Runs the complete bpquant pipeline",
        epilog=EPILOG
    )
    add_pipeline_args(pipeline_parser)

    csv_to_fasta_parser = subparsers.add_parser(
        "csv2fasta",
        description="Converts CSV formatted input sequences to FASTA format",
        epilog=EPILOG,
    )
    add_csv_to_fasta_args(csv_to_fasta_parser)

    indexing_parser = subparsers.add_parser(
        "index",
        description="Generates index for FASTA file",
        epilog=EPILOG,
    )
    add_indexing_args(indexing_parser)

    alignment_parser = subparsers.add_parser(
        "align",
        description="Aligns input reads against custom index",
        epilog=EPILOG,
    )
    add_aligner_args(alignment_parser)

    requantify_parser = subparsers.add_parser(
        "count",
        description="Counts reads aligning to custom reference and classifies them",
        epilog=EPILOG,
    )
    add_requantify_args(requantify_parser)

    plot_reads_parser = subparsers.add_parser(
        "plot",
        description="Plot reads aligning to custom references",
        epilog=EPILOG,
    )
    add_plot_reads_args(plot_reads_parser)

    args = parser.parse_args()

    try:
        args.func(args)
    except AttributeError as e:
        print(e)
        parser.parse_args(["--help"])
