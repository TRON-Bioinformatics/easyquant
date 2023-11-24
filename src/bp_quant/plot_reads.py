from argparse import ArgumentParser
import csv
import sys

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib.path import Path

csv.field_size_limit(sys.maxsize)


def get_seq_to_pos(seq_table_file):
    """
    Parses the sequence table and returns a dict from name to pos
    """

    # seq_to_pos is a dict from sequence name to the position of interest (breakpoint or junction)

    # initialize dict
    seq_to_pos = {}

    with open(seq_table_file, "r") as csvfile:
        # Auto detect dialect of input file
        dialect = csv.Sniffer().sniff(csvfile.readline(), delimiters=";,\t")
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, dialect=dialect)

        # Iterate over input file rows
        for row in reader:
            name = row["name"]

            pos_arr = None
            
            if "," in row["position"]:
                pos_arr = row["position"].split(",")
            else:
                pos_arr = [0, row["position"], len(row["sequence"])]
            intervals = []
            for i in range(len(pos_arr)-1):
                interval_name = "{}_{}".format(pos_arr[i], pos_arr[i+1])
                intervals.append((interval_name, int(pos_arr[i]), int(pos_arr[i+1])))
            seq_to_pos[name] = intervals

    return seq_to_pos


def plot(seq_to_pos, read_info, output_path):
    """

    """
    read_dict = {}

    with open(read_info) as inf:
        next(inf)
        for line in inf:
            elements = line.rstrip().split("\t")
            read_name = elements[0]
            read_group = elements[1]
            seq_name = elements[2]
            start = int(elements[3])
            stop = int(elements[4])
            mismatches = int(elements[5])
            read_type = elements[6]
            
            if seq_name not in read_dict:
                read_dict[seq_name] = {}

            if read_name not in read_dict[seq_name]:
                read_dict[seq_name][read_name] = []

            read_dict[seq_name][read_name].append((read_group, start, stop, read_type))


    pdf = PdfPages(output_path)
        
    #fig, ax = plt.figure(figsize=(3, 3))
        
    for seq_name in seq_to_pos:

        fig, ax = plt.subplots()
        fig.set_size_inches(20, 150)
        h = 10
        min_start = 0
        max_stop = 0
        for interval_name, ref_start, ref_stop in seq_to_pos[seq_name]:
            plt.axvline(x=ref_stop, color='grey', linewidth=1)
            min_start = min(min_start, ref_start)
            max_stop = max(max_stop, ref_stop)
        plt.xticks(range(min_start, max_stop+100, 100))


        #J00128:23:H323KBBXX:6:2117:31084:35194  R1      CLDN18_1        44      95      within
        for read_name in read_dict[seq_name]:
            r1_rg, r1_start, r1_stop, r1_type = read_dict[seq_name][read_name][0]
            r2_rg, r2_start, r2_stop, r2_type = read_dict[seq_name][read_name][1]


            col_dict = {"span": "yellow", "junc": "red", "within": "green", "softjunc": "blue"}
                
            r1 = Rectangle((r1_start, h), r1_stop-r1_start, 4, color=col_dict[r1_type])
            r2 = Rectangle((r2_start, h), r2_stop-r2_start, 4, color=col_dict[r2_type])

            plt.plot([r1_stop+2, r2_start-3], [h+2, h+2], color='grey', linestyle='dashed', linewidth=0.5)

            ax.add_patch(r1)
            ax.add_patch(r2)

            h += 12


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
    seq_to_pos = get_seq_to_pos(args.seq_table_file)

    plot(seq_to_pos, args.input_read_info, args.output)
    
