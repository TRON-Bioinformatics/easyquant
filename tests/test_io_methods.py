"""
Tests for IO module.
"""

import tempfile
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

        with open("test.sam", "w", encoding="utf8") as outf:
            outf.write("@HD\tVN:1.5\tSO:unsorted\tGO:query\t@SQ\tSN:c42c085fa5d881fdd75aad0467853d8e\tLN:800\n")
            outf.write("A00259:290:HKKM5DMXY:1:1266:18267:23735 99      0ef1d60df2e9ba9378f54b4ee3d7da0f        190     1       100M    =       237     147           CTCAGTGGCACACTTCACATTGAAGATGCAGCTGTGCCAGGCGCACAGGAAGCTCTTAAAAGGAATACAAACAAGTGATCCTAATGCTGTGGTCATGGGA    FFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  AS:i:0  XS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0        MD:Z:100        YS:i:0  YT:Z:CP\n")

        self.assertEqual(get_read_count_bam("test.sam"), 1)
