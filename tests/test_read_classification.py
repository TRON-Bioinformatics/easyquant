#!/usr/bin/env python

import unittest

from bp_quant.read_classification import calc_anchor
from bp_quant.read_classification import check_for_indels
from bp_quant.read_classification import classify_read
from bp_quant.read_classification import count_mismatches_in_region
from bp_quant.read_classification import get_match_list
from bp_quant.read_classification import get_query_pos
from bp_quant.read_classification import get_read_type


class TestReadClassification(unittest.TestCase):
    
    def setUp(self):
        self.aln_pairs = [
            (0, 5, 'C'), 
            (1, 6, 't'),
            (2, 7, 'T'),
            (3, 8, 'A'),
            (4, 9, 'A'),
            (5, 10, 'C'),
            (6, 10, 'C'),
            (7, 11, 'T'),
            (8, 12, 'T'),
            (9, 13, 'G'),
            (10, 14, 'A')
        ]


    def test_calc_anchor(self):
        # Positive test case
        res = calc_anchor(aln_start = 88, aln_end = 120, ref_end = 100)
        self.assertEqual(res, 12)

        # Negative test case if aln_start is negative
        with self.assertRaises(ValueError):
            res = calc_anchor(aln_start = -10, aln_end = 20, ref_end = 10)

        # Negative test case if aln_end is negative
        with self.assertRaises(ValueError):
            res = calc_anchor(aln_start = 0, aln_end = -20, ref_end = 10)

        # Negative test case if ref_end is negative
        with self.assertRaises(ValueError):
            res = calc_anchor(aln_start = 0, aln_end = 20, ref_end = -10)

        # Negative test case if aln_start is greater than aln_end
        with self.assertRaises(ValueError):
            res = calc_anchor(aln_start = 30, aln_end = 10, ref_end = 15)
        
        # Negative test case if aln_start is equal to aln_end
        with self.assertRaises(ValueError):
            res = calc_anchor(aln_start = 10, aln_end = 10, ref_end = 15)

        # Negative test case if aln_start and aln_end are smaller than reference stop
        res = calc_anchor(aln_start = 10, aln_end = 11, ref_end = 15)
        self.assertEqual(res, 0)


    def test_get_query_pos(self):

        # Positive test case
        res = get_query_pos(pos = 5, aln_pairs = self.aln_pairs)
        self.assertEqual(res, 0)

        # Positive test case where position has two elements in list
        res = get_query_pos(pos = 10, aln_pairs = self.aln_pairs)
        self.assertEqual(res, 5)

        # Negative test case where searched position is outside of query
        res = get_query_pos(pos = 15, aln_pairs = self.aln_pairs)
        self.assertEqual(res, -1)


    def test_get_match_list(self):

        res = get_match_list(region_start = 6, region_end = 11, aln_pairs = self.aln_pairs)
        self.assertEqual(res, ['t', 'T', 'A', 'A', 'C', 'C'])


    def test_count_mismatches_in_region(self):

        res = count_mismatches_in_region(5, 14, self.aln_pairs)
        self.assertEqual(res, 1)


    def test_check_for_indels(self):
        res = check_for_indels(5, 14, self.aln_pairs)
        self.assertEqual(res, False)


    def test_get_read_type(self):
        # Test case for junction read
        res = get_read_type(10, 12)
        self.assertEqual(res, "junc")

        # Test case for softjunction read
        res = get_read_type(10, 8)
        self.assertEqual(res, "softjunc")

        # Test case for regular read
        res = get_read_type(10, 0)
        self.assertEqual(res, "within")


    def test_classify_read(self):
        intervals = [
            ('0_10', 0, 10),
            ('10_14', 10, 14)
        ]

        # Test case for regular read

        read_info = {
            "class": "within",
            "interval": '10_14',
            "anchor": 0,
            "nm": 0,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": False
        }
        res = classify_read(
            aln_start = 11, 
            aln_end = 14, 
            aln_pairs = self.aln_pairs, 
            intervals = intervals, 
            bp_dist = 3
        )
        self.assertEqual(res, read_info)

        # Test case for softjunct read
        read_info = {
            "class": "softjunc",
            "interval": '0_10',
            "anchor": 4,
            "nm": 1,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": False
        }
        res = classify_read(
            aln_start = 5, 
            aln_end = 14, 
            aln_pairs = self.aln_pairs, 
            intervals = intervals, 
            bp_dist = 5
        )
        self.assertEqual(res, read_info)

        # Test case for junction read with SNPs
        read_info = {
            "class": "junc",
            "interval": '0_10',
            "anchor": 4,
            "nm": 1,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": True
        }

        res = classify_read(
            aln_start = 5, 
            aln_end = 14, 
            aln_pairs = self.aln_pairs, 
            intervals = intervals, 
            bp_dist = 3
        )
        self.assertEqual(res, read_info)

        # Test case for junction read with INDELs
        

        read_info = {
            "class": "junc",
            "interval": '0_10',
            "anchor": 4,
            "nm": 1,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": True
        }

        res = classify_read(
            aln_start = 5, 
            aln_end = 14, 
            aln_pairs = self.aln_pairs, 
            intervals = intervals, 
            bp_dist = 3
        )
        self.assertEqual(res, read_info)