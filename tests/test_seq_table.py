#!/usr/bin/env python

import unittest
from unittest.mock import MagicMock, patch, mock_open

from bp_quant.seq_table import generate_intervals
from bp_quant.seq_table import read_table
from bp_quant.seq_table import get_seq_to_pos_dict


class TestSeqTable(unittest.TestCase):
    
    def setUp(self):
        # Mock seq table file
        # open_mock = mock_open()
        self.SEQ_TABLE_FILE = "test_seq_table.csv"
        # with patch("__main__.open", open_mock, create=True):
        #     with open(self.SEQ_TABLE_FILE, "w") as outf:
        #         outf.write("CLDN18_1;ACGTACGT;0,4,8")


    def tearDown(self):
        # Delete seq table file
        pass

    def test_generate_intervals(self):
        res = generate_intervals([0, 4, 8])
        self.assertEqual(
            res,
            [
                ('0_4', 0, 4),
                ('4_8', 4, 8)
            ]
        )


    def test_read_table(self):
        mock_foo = MagicMock()
        mock_foo.read_table.return_value = read_table(self.SEQ_TABLE_FILE)
        res = list(mock_foo.read_table())
        print(res)
        self.assertEqual(res, ["CLDN18_1;ACGTACGT;0,4,8"])


    def test_get_seq_to_pos(self):
        self.assertEqual(
            get_seq_to_pos_dict(self.SEQ_TABLE_FILE), 
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
