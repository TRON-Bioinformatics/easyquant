#!/usr/bin/env python

import unittest

from bp_quant.coverage import perc_true


class TestRequantify(unittest.TestCase):

    def test_perc_true(self):
        data = [4, 3, 0, 5, 9, 0, 1, 0, 0, 0]
        self.assertEqual(perc_true(data), 0.5)


if __name__ == "__main__":
    unittest.main()
