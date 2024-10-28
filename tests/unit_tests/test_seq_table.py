#!/usr/bin/env python

import tempfile
import unittest

# pylint: disable=E0401
from bp_quant.io.seq_table import generate_intervals
from bp_quant.io.seq_table import read_table
from bp_quant.io.seq_table import get_seq_to_pos_dict


class TestSeqTable(unittest.TestCase):
    """Provides unit tests for seq table module."""

    def setUp(self):
        # Create seq table file
        pass


    def tearDown(self):
        # Delete seq table file
        pass

    def test_generate_intervals(self):
        """Test case for generating intervals."""

        res = generate_intervals([0, 4, 8])
        self.assertEqual(
            res,
            [
                ('0_4', 0, 4),
                ('4_8', 4, 8)
            ]
        )


    @unittest.skip
    def test_read_table(self):
        """Test case for reading the input table."""

        with tempfile.NamedTemporaryFile() as temp:
            temp.write(b"CLDN18_1;ACGTACGT;0,4,8")
            temp.read()
            self.assertEqual(read_table(temp.name), ["CLDN18_1;ACGTACGT;0,4,8"])


    @unittest.skip
    def test_get_seq_to_pos(self):
        """Test case for generating the seq_to_pos dictionary from the input table."""

        with tempfile.NamedTemporaryFile() as temp:
            temp.write(b"CLDN18_1;ACGTACGT;0,4,8")
            temp.read()
            self.assertEqual(
                get_seq_to_pos_dict(temp.name),
                {
                    "CLDN18_1": (
                        [
                            ('0_400', 0, 400),
                            ('400_786', 400, 786)
                        ],
                        "0,400,786"
                    )
                }
            )
