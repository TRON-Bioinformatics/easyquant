"""
IO module to execute commands and get read counts.
"""

import logging
import os
import subprocess
import sys


def create_folder(folder_path: str) -> bool:
    """Creates folder on the filesystem if it doesn't exist."""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        return True
    return False


def execute_cmd(cmd: str, working_dir = ".") -> bool:
    """This function pushes a command into a subprocess."""
    logging.info("Executing CMD: %s", cmd)
    p = subprocess.run(
        cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        cwd = working_dir,
        check=False,
        shell=True
    )
    if p.returncode != 0:
        logging.error("Command \"%s\" returned non-zero exit status", cmd)
        logging.error(p.stderr)
        sys.exit(1)
    return True


def get_read_count(infile: str, input_format: str) -> int:
    """Returns the read count from the input file."""
    if input_format == "fq":
        return get_read_count_fq(infile)
    if input_format == "bam":
        return get_read_count_bam(infile)
    return 0


def get_read_count_fq(fq_file: str) -> int:
    """Parses input FASTQ to get read count"""
    with subprocess.Popen(("zcat", "-f", fq_file), stdout=subprocess.PIPE) as ps:
        result = subprocess.check_output(("wc", "-l"), stdin=ps.stdout)
        return int(int(result) / 4)
    return 0


def get_read_count_bam(bam_file: str) -> int:
    """Parses input BAM to get read count"""
    result = subprocess.check_output(["samtools", "view", "-c", bam_file])
    return int(result)
