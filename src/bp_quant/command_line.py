import argparse
import bp_quant
from bp_quant.csv_to_fasta import add_csv_to_fasta_args
from bp_quant.indexing import add_indexing_args
from bp_quant.aligning import add_aligner_args
from bp_quant.requantify import add_requantify_args
from bp_quant.plot_reads import add_plot_reads_args
from bp_quant.pipeline import add_pipeline_args
from bp_quant.version import version

epilog = "Copyright (c) 2023 TRON gGmbH (See LICENSE for licensing details)"


def bp_quant_cli():
    parser = argparse.ArgumentParser(
        description="bpquant version {}".format(version),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )

    subparsers = parser.add_subparsers(description="Commands")

    pipeline_parser = subparsers.add_parser(
        "pipeline",
        description="Runs the complete bpquant pipeline",
        epilog=epilog
    )
    add_pipeline_args(pipeline_parser)

    csv_to_fasta_parser = subparsers.add_parser(
        "csv2fasta",
        description="Converts CSV formatted input sequences to FASTA format",
        epilog=epilog,
    )
    add_csv_to_fasta_args(csv_to_fasta_parser)

    indexing_parser = subparsers.add_parser(
        "index",
        description="Generates index for FASTA file",
        epilog=epilog,
    )
    add_indexing_args(indexing_parser)

    alignment_parser = subparsers.add_parser(
        "align",
        description="Aligns input reads against custom index",
        epilog=epilog,
    )
    add_aligner_args(alignment_parser)

    requantify_parser = subparsers.add_parser(
        "count",
        description="Counts reads aligning to custom reference and classifies them",
        epilog=epilog,
    )
    add_requantify_args(requantify_parser)

    plot_reads_parser = subparsers.add_parser(
        "plot",
        description="Plot reads aligning to custom references",
        epilog=epilog,
    )
    add_plot_reads_args(plot_reads_parser)

    args = parser.parse_args()

    try:
        args.func(args)
    except AttributeError as e:
        parser.parse_args(["--help"])
