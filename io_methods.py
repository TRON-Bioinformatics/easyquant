import csv
import logging
import os
import subprocess
import sys

csv.field_size_limit(sys.maxsize)


def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

def execute_cmd(cmd, working_dir = "."):
    """This function pushes a command into a subprocess."""
    logging.info("Executing CMD: {}".format(cmd))
    p = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = working_dir, shell=True)
    if p.returncode != 0:
        logging.error("Command \"{}\" returned non-zero exit status".format(cmd))
        logging.error(p.stderr)
        sys.exit(1)


def csv_to_fasta(csv_in, fasta_out):
    """This function converts the target sequences TSV/CSV file to the FASTA format."""
    outf = open(fasta_out, "w")
    with open(csv_in, "r", newline="\n") as csvfile:
        # Auto detect dialect of input file
        dialect = csv.Sniffer().sniff(csvfile.readline(), delimiters=";,\t")
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, dialect=dialect)

        # Iterate over input file rows
        for row in reader:

            name = row["name"]
            sequence = row["sequence"]
            position = row["position"]

            outf.write(">{}\n".format(name))
            for i in range(0, len(sequence), 60):
                outf.write("{}\n".format(sequence[i:i+60]))

    outf.close()


def get_fasta_size(fasta):
    """This function calculates the FASTA size in bp."""
    fasta_size = 0
    with open(fasta) as inf:
        for line in inf:
            if not line.rstrip().startswith(">"):
                fasta_size += len(line.rstrip())
    return fasta_size
