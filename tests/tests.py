#!/usr/bin/env python

import os
import unittest
import pysam

SEQ_TABLE_FILE = os.path.join("example_data", "CLDN18_Context_seq.csv")
BAM_FILE = os.path.join("example_data", "example_rna-seq.bam")

from bp_quant.requantify import Quantification, process_secondary_alignments, get_aligner, perc_true, get_seq_to_pos, classify_read, is_chimeric_alignment, is_singleton


class TestRequantify(unittest.TestCase):

    def test_perc_true(self):
        data = [4, 3, 0, 5, 9, 0, 1, 0, 0, 0]
        self.assertEqual(perc_true(data), 0.5)


    def test_process_secondary_alignments(self):
        read_dict = {
            'ref_name_1': {
                'query_name_1': {
                    'R1': [{
                        'query_name': 'query_name_1', 
                        'first_in_pair': True, 
                        'unmapped': False, 
                        'reference_name': 'ref_name_1', 
                        'flag': 355, 
                        'cigar': '100M', 
                        'start': 141, 
                        'stop': 241, 
                        'pairs': []
                    }], 
                    'R2': [{
                        'query_name': 'query_name_1', 
                        'first_in_pair': False, 
                        'unmapped': False, 
                        'reference_name': 'ref_name_1', 
                        'flag': 403, 
                        'cigar': '100M', 
                        'start': 233, 
                        'stop': 333, 
                        'pairs': []
                    }]
                }
            }, 
            'ref_name_2': {
                'query_name_2': {
                    'R1': [], 
                    'R2': [{
                        'query_name': 'query_name_2', 
                        'first_in_pair': False, 
                        'unmapped': False, 
                        'reference_name': 'ref_name_2', 
                        'flag': 401, 
                        'cigar': '100M', 
                        'start': 284, 
                        'stop': 384, 
                        'pairs': []
                    }]
                }
            }, 
            'ref_name_3': {
                'query_name_2': {
                    'R1': [{
                        'query_name': 'query_name_2', 
                        'first_in_pair': True, 
                        'unmapped': False, 
                        'reference_name': 'ref_name_3', 
                        'flag': 353, 
                        'cigar': '100M', 
                        'start': 23, 
                        'stop': 123, 
                        'pairs': []
                    }], 
                    'R2': []
                }
            }
        }
        read_pairings = [
            ({
                'query_name': 'query_name_1', 
                'first_in_pair': True, 
                'unmapped': False, 
                'reference_name': 'ref_name_1', 
                'flag': 355, 
                'cigar': '100M', 
                'start': 141, 
                'stop': 241, 
                'pairs': []
            },
            {
                'query_name': 'query_name_1',
                'first_in_pair': False,
                'unmapped': False,
                'reference_name': 'ref_name_1',
                'flag': 403,
                'cigar': '100M',
                'start': 233,
                'stop': 333,
                'pairs': []
            }),
            ({
                "reference_name": "ref_name_2",
                "query_name": "query_name_2",
                "unmapped": True,
                "flag": 325,
                "start": -1,
                "stop": -1,
                "pairs": None,
                "cigar": None
            },
            {
                'query_name': 'query_name_2', 
                'first_in_pair': False, 
                'unmapped': False, 
                'reference_name': 'ref_name_2', 
                'flag': 401, 
                'cigar': '100M', 
                'start': 284, 
                'stop': 384, 
                'pairs': []
            }),
            ({
                'query_name': 'query_name_2', 
                'first_in_pair': True, 
                'unmapped': False, 
                'reference_name': 'ref_name_3', 
                'flag': 353, 
                'cigar': '100M', 
                'start': 23, 
                'stop': 123, 
                'pairs': []
            },
            {
                "reference_name": "ref_name_3",
                "query_name": "query_name_2",
                "unmapped": True,
                "flag": 389,
                "start": -1,
                "stop": -1,
                "pairs": None,
                "cigar": None
            })
        ]
        self.assertEqual(process_secondary_alignments(read_dict), read_pairings)


if __name__ == "__main__":
    unittest.main()
