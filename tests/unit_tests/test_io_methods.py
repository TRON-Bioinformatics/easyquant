"""
Tests for IO module.
"""

import tempfile
import unittest

# pylint: disable=E0401
from bp_quant.io.io_methods import create_folder
from bp_quant.io.io_methods import execute_cmd
from bp_quant.io.io_methods import get_read_count_fq
from bp_quant.io.io_methods import get_read_count_bam

class TestIOMethods(unittest.TestCase):
    """Provides unit tests for IO module."""

    @unittest.skip
    def test_create_folder(self):
        """
        Test case for folder creation.
        """

        self.assertFalse(create_folder("test"))


    def test_execute_command(self):
        """
        Test case for command execution.
        """
        self.assertTrue(execute_cmd("ls"))


    def test_get_read_count_fq(self):
        """
        Test case for read counting in a FASTQ file.
        """
        with tempfile.NamedTemporaryFile() as temp:
            temp.write(b"@Read1\nACGT\n+\n;;;;\n")
            temp.read()
            self.assertEqual(get_read_count_fq(temp.name), 1)


    @unittest.skip
    def test_get_read_count_bam(self):
        """
        Test case for read counting in a BAM file.
        """

        with tempfile.NamedTemporaryFile() as temp:
            temp.write("@HD\n")
            temp.write("test_read 99 XY 190 1 100M = 237 147 ACGT FFFF\n")
            temp.read()
            self.assertEqual(get_read_count_bam("test.sam"), 1)
