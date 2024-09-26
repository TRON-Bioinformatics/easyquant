"""
Tests for IO module.
"""


import unittest

from bp_quant.io_methods import create_folder
from bp_quant.io_methods import execute_cmd
from bp_quant.io_methods import get_read_count_fq
from bp_quant.io_methods import get_read_count_bam

class TestReadProcessing(unittest.TestCase):
    """Provides unit tests for IO module."""

    def test_create_folder(self):
        """
        Test case for folder creation.
        """
        self.assertTrue(create_folder("test"))


    def test_execute_command(self):
        """
        Test case for command execution.
        """
        self.assertTrue(execute_cmd("ls"))


    def test_get_read_count_fq(self):
        """
        Test case for read counting in a FASTQ file.
        """
        with open("test.fq", "w", encoding="utf8") as outf:
            outf.write("@Read1\nACGT\n+\n;;;;")

        self.assertEqual(get_read_count_fq("test.fq"), 1)


    def test_get_read_count_bam(self):
        """
        Test case for read counting in a BAM file.
        """

        with open("test.sam", "w", encoding="utf8") as outf:
            outf.write("Mocking an alignment entry.")

        self.assertEqual(get_read_count_bam("test.sam"), 1)
