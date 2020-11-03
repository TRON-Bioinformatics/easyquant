#!/usr/bin/env python

"""
Junction/Spanning read pair counting

@author: BNT (URLA)
@version: 20190118
"""

from __future__ import print_function, division
import os.path
from argparse import ArgumentParser
#from urla_logger_latest import Logger
import pysam # pysam is not available for windows (where I run pylint) => pylint: disable=E0401


# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class Requantification(object):
    """Select alignments belonging to putative fusions from an s/bam file"""
    def __init__(self, bam, bed, output, bp_distance_threshold):
        """Parameter initialization"""
        self.in_bam = pysam.AlignmentFile(bam, "rb")
        self.bed = bed
        self.output = output
        self.bp_distance_threshold = int(bp_distance_threshold)
        self.fusion_seq_dict = self.get_fusions()
        #self.logger = Logger("{}.fusionReadFilterLog".format(output))
        self.input_read_count = 0

    def get_fusions(self):
        fusion_seq_dict = {}
        with open(self.bed) as bed:
            for line in bed:
                elements = line.rstrip().split("\t")
                fgid = elements[0]
                bp = elements[1]

                fusion_seq_dict[fgid] = [int(bp), 0, 0, 0, 0, 0]

        return fusion_seq_dict

    def quantify_read_groups(self, read_buffer, reference_base):
        """Count junction/spanning pairs on ft and wt1/2"""
        anchor = 0
        # the read_buffer contains all reads mapping to the fusion transcript and the wildtyp backgrounds
        # the read_group is list of reads (incl multimappers) belonging to a single read pair
        # reads from this read group can map independently to ft, wt1 and wt2
        for read_group in read_buffer:
            junction = 0 # 0 = no junction overlapping read, 1+ = at least one of the paired reads is a junction read
            spanning = [0, 0] # left: paired read cnt fully left of the fusion bp, right: paired read cnt fulls right of the fusion bp

            # each processed read belongs to one of the three defined read_groups
            for read in read_buffer[read_group]:
                junction, spanning, anchor = self.count_junc_span(read, self.fusion_seq_dict[reference_base][0], junction, spanning, anchor)

            #print("stored before: {}".format(self.fusion_seq_dict[reference_base]))
            #print("jft: {}, sft: {}, jwt: {}, swt: {}".format(is_junction_ft, is_spanning_ft, is_junction_wt, is_spanning_wt))
            # update junction and spanning counts
            self.update_counts(reference_base, junction, spanning, 1)
        # set max anchor for this reference base
        self.fusion_seq_dict[reference_base][5] = anchor
            #print("stored after: {}".format(self.fusion_seq_dict[reference_base]))

    def count_junc_span(self, read, breakpoint_pos, junction_count, spanning_counts, anchor):
        """Identify whether a read belongs to a junction or spanning pair"""
        # read overlaps the fusion breakpoint
        if breakpoint_pos >= self.bp_distance_threshold and read.get_overlap(breakpoint_pos - self.bp_distance_threshold, breakpoint_pos + self.bp_distance_threshold) == 2 * self.bp_distance_threshold:
            junction_count += 1
            anchor = max(anchor, min(read.reference_end - breakpoint_pos, breakpoint_pos - read.reference_start))
        # read maps left of the fusion breakpoint
        if read.reference_end < (breakpoint_pos + self.bp_distance_threshold):
            spanning_counts[0] += 1
        # read maps right of the fusion breakpoint
        if read.reference_start > (breakpoint_pos - self.bp_distance_threshold):
            spanning_counts[1] += 1
        return junction_count, spanning_counts, anchor

    # counter_start should be 1 for ft, 6 for wt1 and 11 for wt2 (see comment in header_to_dict for details)
    def update_counts(self, reference_base, junctions, spannings, counter_start):
        """Updates the read counter in the fusion_seq_dict"""
        if spannings[0] > 0:
            self.fusion_seq_dict[reference_base][1] += 1
        if spannings[1] > 0:
            self.fusion_seq_dict[reference_base][2] += 1
        if junctions > 0:
            self.fusion_seq_dict[reference_base][3] += 1
        elif spannings[0] > 0 and spannings[0] == spannings[1]:
            self.fusion_seq_dict[reference_base][4] += 1

    def normalize_counts_cpm(self, count):
        """Returns the CPM of a read count based on the number of original input reads"""
        if self.input_read_count == 0:
            return count
        return count / self.input_read_count * 1000000

    def run(self):
        """Walk linewise through a s/bam file and send reads mapping to the same fusion/wt context to counting"""
        #self.logger.info("Starting fusion read filtering")
        count_lines = 0
        count_processed_refs = 0
        # the read buffer is a dict of lists. The dict has query names as keys and a list of corresponding reads as values.
        read_buffer = {}
        last_reference = ""

        for read in self.in_bam.fetch():
            count_lines += 1
            if not last_reference in read.reference_name and count_lines > 1:
                self.quantify_read_groups(read_buffer, last_reference)
                count_processed_refs += 1
                if count_processed_refs % (len(self.fusion_seq_dict) / 10) == 0:
                    print("At least {}% completed...".format(count_processed_refs / (len(self.fusion_seq_dict) / 100)))
                read_buffer.clear()
            if not read.query_name in read_buffer:
                read_buffer[read.query_name] = []
            read_buffer[read.query_name].append(read)
            last_reference = read.reference_name


        # process the last read_buffer
        self.quantify_read_groups(read_buffer, last_reference)
        print("All done! Reads did not map to ~{}% of the reference input".format(100 - (count_processed_refs / (len(self.fusion_seq_dict) / 100))))
        # close reading/writing stream
        self.in_bam.close()
        # write data to file. Change out_file_sep to ";" to create a csv file

        # get input read count
        with open(os.path.join(os.path.dirname(self.output), "star", "Log.final.out"), "r") as rfile:
            for i, line in enumerate(rfile):
                elements = line.rstrip().split()
                if i == 5:
                    self.input_read_count = int(elements[-1])

        out_file_sep = ";"
        header_string = out_file_sep.join(
            ["name", "position", "a", "b", "junc", "span", "anch"])
        # write counts and normalized counts
        # for normalization, breakpoint and anchor must not be converted: list positions 0, 5, 10, 15
        no_norm = {0, 5, 6, 11, 12, 17}
        with open("{}.counts".format(self.output), "w") as out_counts, open(self.output, "w") as out_norm:
            out_counts.write("{}\n".format(header_string))
            out_norm.write("{}\n".format(header_string))
            for key in self.fusion_seq_dict:
                out_counts.write("{}{}".format(key, out_file_sep))
                out_norm.write("{}{}".format(key, out_file_sep))
                # write counts directly
                out_counts.write("{}\n".format(out_file_sep.join(map(str, self.fusion_seq_dict[key]))))
                # normalize and write normalized counts
                self.fusion_seq_dict[key] = [read_count if i in no_norm else self.normalize_counts_cpm(read_count) for i, read_count in enumerate(self.fusion_seq_dict[key])]
                #self.fusion_seq_dict[key] = map(self.normalize_counts_cpm, self.fusion_seq_dict[key])
                out_norm.write("{}\n".format(out_file_sep.join(map(str, self.fusion_seq_dict[key]))))

def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input', dest='input', help='Specify input BAM file', required=True)
    parser.add_argument('-b', '--bed', dest='bed', help='Specify input BED file', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Specify output file', required=True)
    parser.add_argument('-d', '--bp_distance', dest='bp_distance', help='Threshold of bases around the breakpoint for junction/spanning counting')
    args = parser.parse_args()

    requant = Requantification(args.input, args.bed, args.output, args.bp_distance)
    requant.run()

if __name__ == '__main__':
    main()
