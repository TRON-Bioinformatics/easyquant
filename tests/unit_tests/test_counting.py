"""
Tests for read classification module.
"""

import unittest

# pylint: disable=E0401
from bp_quant.quantification.counting import init_count_dict_interval_mode
from bp_quant.quantification.counting import init_count_dict_regular_mode
from bp_quant.quantification.counting import init_cov_dict
from bp_quant.quantification.counting import get_spanning_intervals
from bp_quant.quantification.counting import count_reads_interval_mode
from bp_quant.quantification.counting import count_reads_regular_mode


class TestCounting(unittest.TestCase):
    """Provides unit tests for read classification module."""

    def setUp(self):
        self.seq_to_pos_dict = {
            "CLDN18_1": (
                [
                    ('0_400', 0, 400),
                    ('400_786', 400, 786)
                ],
                "0,400,786"
            )
        }
        self.intervals = [
            ('0_400', 0, 400),
            ('400_786', 400, 786)
        ]


    def test_init_count_dict_interval_mode(self):
        """Test case for count_dict_interval_mode."""

        res = init_count_dict_interval_mode(self.seq_to_pos_dict)
        exp = {
            'CLDN18_1': {
                '0_400': {
                    'coverage_mean': 0.0,
                    'coverage_median': 0.0,
                    'coverage_perc': 0.0,
                    'overlap_interval_end_reads': 0,
                    'span_interval_end_pairs': 0,
                    'within_interval': 0
                },
                '400_786': {
                    'coverage_mean': 0.0,
                    'coverage_median': 0.0,
                    'coverage_perc': 0.0,
                    'overlap_interval_end_reads': 0,
                    'span_interval_end_pairs': 0,
                    'within_interval': 0
                }
            }
        }
        self.assertEqual(res, exp)


    def test_init_count_dict_regular_mode(self):
        """Test case for count_dict_regular_mode."""
        res = init_count_dict_regular_mode(self.seq_to_pos_dict)
        exp = {
            'CLDN18_1': {
                'junc': 0,
                'span': 0,
                'anch': 0,
                'a': 0,
                'b': 0
            }
        }
        self.assertEqual(res, exp)


    @unittest.skip
    def test_init_cov_dict(self):
        """Test case for init_cov_dict."""
        res = init_cov_dict(self.seq_to_pos_dict)
        exp = {}
        self.assertEqual(res, exp)


    def test_get_spanning_intervals(self):
        """Test case for get_spanning_intervals."""
        intervals = [
            ('0_400', 0, 400),
            ('400_786', 400, 786),
            ('786_800', 786, 800),
            ('800_1200', 800, 1200)
        ]

        res = get_spanning_intervals(intervals, "0_400", "786_800")
        exp = ['0_400', '400_786']
        self.assertEqual(res, exp)

        res = get_spanning_intervals(intervals, "400_786", "800_1200")
        exp = ['400_786', '786_800']
        self.assertEqual(res, exp)


    def test_count_reads_interval_mode(self):
        """Test case for count_reads_interval_mode."""
        r1_info = {
            "class": "within",
            "interval": '400_786',
            "anchor": 0,
            "nm": 0,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": False
        }
        r2_info = {
            "class": "junc",
            "interval": '0_400',
            "anchor": 38,
            "nm": 2,
            "nm_in_bp_area": 2,
            "contains_snp_or_indel_in_bp_area": True
        }
        res = count_reads_interval_mode(r1_info, r2_info, self.intervals, True)
        exp = {
            '0_400': {
                'overlap_interval_end_reads': 1,
                'span_interval_end_pairs': 0,
                'within_interval': 0
            },
            '400_786': {
                'overlap_interval_end_reads': 0,
                'span_interval_end_pairs': 0,
                'within_interval': 1
            }
        }
        self.assertEqual(res, exp)


    def test_count_reads_regular_mode(self):
        """Test case for count_reads_regular_mode."""

        # Junction test case
        r1_info = {
            "class": "within",
            "interval": '400_786',
            "anchor": 0,
            "nm": 0,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": False
        }
        r2_info = {
            "class": "junc",
            "interval": '0_400',
            "anchor": 38,
            "nm": 2,
            "nm_in_bp_area": 2,
            "contains_snp_or_indel_in_bp_area": True
        }
        res = count_reads_regular_mode(r1_info, r2_info, self.intervals, True)
        exp = {'junc': 1, 'span': 0, 'a': 0, 'b': 1}
        self.assertEqual(res, exp)

        # Spanning pair test case
        r1_info = {
            "class": "span",
            "interval": '400_786',
            "anchor": 0,
            "nm": 0,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": False
        }
        r2_info = {
            "class": "span",
            "interval": '0_400',
            "anchor": 0,
            "nm": 2,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": True
        }
        res = count_reads_regular_mode(r1_info, r2_info, self.intervals, True)
        exp = {'junc': 0, 'span': 1, 'a': 1, 'b': 1}
        self.assertEqual(res, exp)