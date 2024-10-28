"""
This module generates a plot from the read info table.
"""

# pylint: disable=E0401
from matplotlib.backends.backend_pdf import PdfPages # type: ignore
import matplotlib.pyplot as plt # type: ignore
from matplotlib.patches import Rectangle # type: ignore

from bp_quant.io.seq_table import get_seq_to_pos_dict


def parse_read_info(read_info_file: str) -> dict:
    """
    Parses the read info file and stores the information into a dictionary.
    """
    read_dict = {}

    with open(read_info_file, encoding="utf8") as inf:
        next(inf)
        for line in inf:
            elements = line.rstrip().split("\t")
            read_name = elements[0]
            read_group = elements[1]
            seq_name = elements[2]
            start = int(elements[3])
            stop = int(elements[4])
            read_type = elements[6]

            if seq_name not in read_dict:
                read_dict[seq_name] = {}

            if read_name not in read_dict[seq_name]:
                read_dict[seq_name][read_name] = []

            read_dict[seq_name][read_name].append((read_group, start, stop, read_type))
    return read_dict


def plot_reads(read_dict: dict, seq_name: str, ax: object) -> None:
    """
    Plots the individual reads.
    """
    h = 10
    for read_name in read_dict[seq_name]:
        _, r1_start, r1_stop, r1_type = read_dict[seq_name][read_name][0]
        _, r2_start, r2_stop, r2_type = read_dict[seq_name][read_name][1]


        col_dict = {"span": "yellow", "junc": "red", "within": "green", "softjunc": "blue"}

        r1 = Rectangle((r1_start, h), r1_stop-r1_start, 4, color=col_dict[r1_type])
        r2 = Rectangle((r2_start, h), r2_stop-r2_start, 4, color=col_dict[r2_type])

        plt.plot(
            [r1_stop+2, r2_start-3],
            [h+2, h+2],
            color='grey',
            linestyle='dashed',
            linewidth=0.5
        )

        ax.add_patch(r1)
        ax.add_patch(r2)

        h += 12

def plot(seq_to_pos: dict, read_info_file: str, output_path: str) -> None:
    """
    This method iterates over the read info file and plots the reads incrementally.
    """

    read_dict = parse_read_info(read_info_file)

    pdf = PdfPages(output_path)

    #fig, ax = plt.figure(figsize=(3, 3))

    for seq_name in seq_to_pos:

        fig, ax = plt.subplots()
        fig.set_size_inches(20, 150)

        min_start = 0
        max_stop = 0
        for _, ref_start, ref_stop in seq_to_pos[seq_name]:
            plt.axvline(x=ref_stop, color='grey', linewidth=1)
            min_start = min(min_start, ref_start)
            max_stop = max(max_stop, ref_stop)
        plt.xticks(range(min_start, max_stop+100, 100))

        plot_reads(read_dict, seq_name, ax)

        plt.title(seq_name)
        pdf.savefig(fig, dpi=500)
        plt.close()
    pdf.close()


def add_plot_reads_args(parser):
    """Add arguments for read plotting to a parser"""
    parser.add_argument(
        "-i",
        "--input_read_info",
        dest="input_read_info",
        help="Input read info file",
        required=True
    )
    parser.add_argument(
        "-t",
        "--seq_table",
        dest="seq_table_file",
        help="Path to input sequence table",
        required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Output file",
        required=True
    )
    parser.set_defaults(func=plot_reads_command)


def plot_reads_command(args):
    """Run read plotting from command line"""
    # parse sequence names and position of interest
    seq_to_pos = get_seq_to_pos_dict(args.seq_table_file)

    plot(seq_to_pos, args.input_read_info, args.output)
