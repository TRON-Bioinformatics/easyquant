"""
File handler module.
"""

# pylint: disable=E0401
import bp_quant.io.file_headers as headers

def write_line_to_file(file_handle: object, values: list) -> bool:
    """Writes a specific line to an open file handle."""
    try:
        file_handle.write(
            "{}\n".format(
                "\t".join([str(e) for e in values])
            ).encode()
        )
        return True
    except IOError:
        return False


def save_plot_script(plot_reads_script: str, read_info_file: str,
                     seq_tab_file: str, plotted_reads_pdf: str,
                     plot_script: str):
    """Saves the plotting command to a shell script.

    Args:
        plot_reads_script (str): Path to the python script for plotting
        read_info_file (str): Path to the read info table
        seq_tab_file (str): Path to the reference table
        plotted_reads_pdf (str): Path to the plot in PDF format
        plot_script (str): Path to the shell script
    """

    # Create plotting script
    cmd = f"{plot_reads_script} \
        -i {read_info_file} \
        -t {seq_tab_file} \
        -o {plotted_reads_pdf}"
    with open(plot_script, "w", encoding="utf8") as outf:
        outf.write(f"#/bin/bash\n\n{cmd}")


def save_counts(output_file: str, seq_to_pos: dict, counts: dict, interval_mode: bool):
    """Saves the count values to a file.

    Args:
        output_file (str): Output file to save the counts to
        seq_to_pos (dict): Interval information
        counts (dict): Count information
        interval_mode (bool): Select if interval mode was used
    """

    with open(output_file, "w", encoding="utf8") as out_handle:
        # write header line
        sp_out = None
        if interval_mode:
            sp_out = headers.COUNTS_INT_MODE
        else:
            sp_out = headers.COUNTS_SINGLE_MODE
        out_line = "\t".join(sp_out) + "\n"
        out_handle.write(out_line)

        # Iterate over sequence dictionary
        for name, _ in seq_to_pos.items():
            # get sequence counts
            seq_counts = counts[name]

            if interval_mode:
                for interval_name in seq_counts:
                    # drop sequence and construct output columns
                    sp_out = [name, interval_name]

                    for colname in headers.COUNTS_INT_MODE[2:]:
                        sp_out.append(str(seq_counts[interval_name][colname]))

                    # write as otput line to output file
                    out_line = "\t".join(sp_out) + "\n"
                    out_handle.write(out_line)
            else:
                position = str(seq_to_pos[name][0][0][2])
                # drop sequence and construct output columns
                sp_out = [name, position]
                for colname in headers.COUNTS_SINGLE_MODE[2:]:
                    sp_out.append(str(seq_counts[colname]))

                # write as otput line to output file
                out_line = "\t".join(sp_out) + "\n"
                out_handle.write(out_line)
