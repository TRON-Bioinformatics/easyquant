#!/usr/bin/env python

"""
This module contains tests for alignment info gathering.
"""

import os
import unittest

# pylint: disable=E0401
import pysam # type: ignore

from bp_quant.validation.alignment_info import get_aligner
from bp_quant.validation.alignment_info import get_sorting
from bp_quant.validation.alignment_info import is_chimeric_alignment
from bp_quant.validation.alignment_info import is_singleton

BAM_FILE = os.path.join("example_data", "example_rna-seq.bam")

class TestAlignmentInfo(unittest.TestCase):
    """Provides unit tests for IO module."""

    def setUp(self):
        # Set up alignment object
        self.aln_obj = pysam.AlignmentFile(BAM_FILE, "rb")


    def test_get_aligner(self):
        """Test get_aligner."""
        self.assertEqual(get_aligner(self.aln_obj), "star")


    def test_get_sorting(self):
        """Test get_sorting."""
        self.assertEqual(get_sorting(self.aln_obj), "unsorted")


    def test_is_chimeric(self):
        """Test is_chimeric."""
        # Test that chimeric read is correctly identified
        read = pysam.AlignedSegment()
        read.query_name = "read_28833_29006_6945"
        read.flag = 99
        read.reference_id = 0
        read.reference_start = 32
        read.next_reference_id = 3
        self.assertTrue(is_chimeric_alignment(read))

        # Check that aligned read with unmapped mapped is identified as valid alignment
        read = pysam.AlignedSegment()
        read.query_name = "read_28833_29006_6945"
        read.flag = 105
        read.reference_id = 0
        read.reference_start = 32
        read.next_reference_id = -1
        self.assertFalse(is_chimeric_alignment(read))

        # Check that unaligned read with mapped mate is identified as valid alignment
        read = pysam.AlignedSegment()
        read.query_name = "read_28833_29006_6945"
        read.flag = 101
        read.reference_id = -1
        read.reference_start = 32
        read.next_reference_id = 3
        self.assertFalse(is_chimeric_alignment(read))


    def test_is_singleton(self):
        """Test is_singleton."""

        # Check that aligned read with unmapped mapped is identified as valid alignment
        read = pysam.AlignedSegment()
        read.query_name = "read_28833_29006_6945"
        read.flag = 105
        read.reference_id = 0
        read.reference_start = 32
        read.next_reference_id = -1
        self.assertTrue(is_singleton(read))

        # Check that unaligned read with mapped mate is identified as valid alignment
        read = pysam.AlignedSegment()
        read.query_name = "read_28833_29006_6945"
        read.flag = 101
        read.reference_id = -1
        read.reference_start = 32
        read.next_reference_id = 3
        self.assertTrue(is_singleton(read))
