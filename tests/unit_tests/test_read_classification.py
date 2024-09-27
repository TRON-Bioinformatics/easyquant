"""
Tests for read classification module.
"""

import unittest

# pylint: disable=E0401
from bp_quant.quantification.read_classification import calc_anchor
from bp_quant.quantification.read_classification import check_for_indels
from bp_quant.quantification.read_classification import classify_read
from bp_quant.quantification.read_classification import count_mismatches_in_region
from bp_quant.quantification.read_classification import get_match_list
from bp_quant.quantification.read_classification import get_query_pos
from bp_quant.quantification.read_classification import get_read_type


class TestReadClassification(unittest.TestCase):
    """Provides unit tests for read classification module."""
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
        """Test calc_anchor."""
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
        """Test get_query_pos."""
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
        """Test get_match_list."""
        res = get_match_list(region_start = 6, region_end = 11, aln_pairs = self.aln_pairs)
        self.assertEqual(res, ['t', 'T', 'A', 'A', 'C', 'C'])

        res = get_match_list(region_start = 6, region_end = 11, aln_pairs = None)
        self.assertEqual(res, [])


    def test_count_mismatches_in_region(self):
        """Test count_mismatches_in_region."""
        res = count_mismatches_in_region(5, 14, self.aln_pairs)
        self.assertEqual(res, 1)


    def test_check_for_indels(self):
        """Test check_for_indels."""
        res = check_for_indels(5, 14, self.aln_pairs)
        self.assertEqual(res, False)


    def test_get_read_type(self):
        """Test get_read_type."""
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
        """Test classify_read."""
        intervals = [
            ('0_10', 0, 10),
            ('10_14', 10, 14)
        ]

        # Test case for regular read

        aln_pairs = [
            (7, 11, 'T'),
            (8, 12, 'T'),
            (9, 13, 'G')
        ]

        read_info = {
            "class": "within",
            "interval": '10_14',
            "anchor": 0,
            "nm": 0,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": False
        }
        res = classify_read(
            aln_pairs = aln_pairs,
            intervals = intervals,
            bp_dist = 3
        )
        self.assertEqual(res, read_info)

        # Test case for softjunc read
        aln_pairs = [
            (0, 5, 'C'),
            (1, 6, 't'),
            (2, 7, 'T'),
            (3, 8, 'A'),
            (4, 9, 'A'),
            (5, 10, 'C'),
            (6, 10, 'C'),
            (7, 11, 'T'),
            (8, 12, 'T'),
            (9, 13, 'G')
        ]

        read_info = {
            "class": "softjunc",
            "interval": '0_10',
            "anchor": 4,
            "nm": 1,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": False
        }
        res = classify_read(
            aln_pairs = aln_pairs,
            intervals = intervals,
            bp_dist = 5
        )
        self.assertEqual(res, read_info)

        # Test case for junction read with SNPs

        aln_pairs = [
            (0, 5, 'C'),
            (1, 6, 't'),
            (2, 7, 'T'),
            (3, 8, 'A'),
            (4, 9, 'A'),
            (5, 10, 'C'),
            (6, 10, 'C'),
            (7, 11, 'T'),
            (8, 12, 'T'),
            (9, 13, 'G')
        ]

        read_info = {
            "class": "junc",
            "interval": '0_10',
            "anchor": 4,
            "nm": 1,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": True
        }

        res = classify_read(
            aln_pairs = aln_pairs,
            intervals = intervals,
            bp_dist = 3
        )
        self.assertEqual(res, read_info)

        # Test case for junction read with INDELs


        read_info = {
            "class": "junc",
            "interval": '0_10',
            "anchor": 5,
            "nm": 1,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": True
        }

        res = classify_read(
            aln_pairs = self.aln_pairs,
            intervals = intervals,
            bp_dist = 3
        )
        self.assertEqual(res, read_info)


        intervals = [
            ('0_272', 0, 272),
            ('272_357', 272, 357)
        ]


        aln_pairs = [
            (0, 210, 'C'),
            (1, 211, 'T'),
            (2, 212, 'C'),
            (3, 213, 'A'),
            (4, 214, 'G'),
            (5, 215, 'T'),
            (6, 216, 'G'),
            (7, 217, 'G'),
            (8, 218, 'C'),
            (9, 219, 'A'),
            (10, 220, 'C'),
            (11, 221, 'A'),
            (12, 222, 'C'),
            (13, 223, 'T'),
            (14, 224, 'T'),
            (15, 225, 'C'),
            (16, 226, 'A'),
            (17, 227, 'C'),
            (18, 228, 'A'),
            (19, 229, 'T'),
            (20, 230, 'T'),
            (21, 231, 'G'),
            (22, 232, 'A'),
            (23, 233, 'A'),
            (24, 234, 'G'),
            (25, 235, 'A'),
            (26, 236, 'T'),
            (27, 237, 'G'),
            (28, 238, 'C'),
            (29, 239, 'A'),
            (30, 240, 'G'),
            (31, 241, 'C'),
            (32, 242, 'T'),
            (33, 243, 'G'),
            (34, 244, 'T'),
            (35, 245, 'G'),
            (36, 246, 'C'),
            (37, 247, 'C'),
            (38, 248, 'A'),
            (39, 249, 'G'),
            (40, 250, 'G'),
            (41, 251, 'C'),
            (42, 252, 'G'),
            (43, 253, 'C'),
            (44, 254, 'A'),
            (45, 255, 'C'),
            (46, 256, 'A'),
            (47, 257, 'G'),
            (48, 258, 'G'),
            (49, 259, 'A'),
            (50, 260, 'A'),
            (51, 261, 'G'),
            (52, 262, 'C'),
            (53, 263, 'T'),
            (54, 264, 'C'),
            (55, 265, 'T'),
            (56, 266, 'T'),
            (57, 267, 'A'),
            (58, 268, 'A'),
            (59, 269, 'A'),
            (60, 270, 'A'),
            (61, 271, 'G'),
            (62, 272, 'G'),
            (63, 273, 'A'),
            (64, 274, 'A'),
            (65, 275, 'T'),
            (66, 276, 'A'),
            (67, 277, 'C'),
            (68, 278, 'A'),
            (69, 279, 'A'),
            (70, 280, 'A'),
            (71, 281, 'C'),
            (72, 282, 'A'),
            (73, 283, 'A'),
            (74, 284, 'G'),
            (75, 285, 'T'),
            (76, 286, 'G'),
            (77, 287, 'A'),
            (78, 288, 'T'),
            (79, 289, 'C'),
            (80, 290, 'C'),
            (81, 291, 'T'),
            (82, 292, 'A'),
            (83, 293, 'A'),
            (84, 294, 'T'),
            (85, 295, 'G'),
            (86, 296, 'C'),
            (87, 297, 'T'),
            (88, 298, 'G'),
            (89, 299, 'T'),
            (90, 300, 'G'),
            (91, 301, 'G'),
            (92, 302, 'T'),
            (93, 303, 'C'),
            (94, 304, 'A'),
            (95, 305, 'T'),
            (96, 306, 'G'),
            (97, 307, 'G'),
            (98, 308, 'G'),
            (99, 309, 'A')
        ]

        read_info = {
            "class": "junc",
            "interval": '0_272',
            "anchor": 38,
            "nm": 0,
            "nm_in_bp_area": 0,
            "contains_snp_or_indel_in_bp_area": False
        }

        res = classify_read(
            aln_pairs = aln_pairs,
            intervals = intervals,
            bp_dist = 10
        )
        self.assertEqual(res, read_info)
