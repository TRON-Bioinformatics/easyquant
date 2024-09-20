"""
Tests for read processing module.
"""


import unittest

from bp_quant.read_processing import process_secondary_alignments

class TestReadProcessing(unittest.TestCase):
    """Provides unit tests for read processing module."""

    def setUp(self):
        self.read_dict = {
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


    def test_process_secondary_alignments(self):
        """Tests for secondary alignment processing."""
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
        self.assertEqual(process_secondary_alignments(self.read_dict), read_pairings)
