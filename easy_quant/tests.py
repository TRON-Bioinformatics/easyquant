#!/usr/bin/env python

import unittest

from easy_quant.requantify import perc_true, mean, median


class TestRequantify(unittest.TestCase):

    def test_perc_true(self):
        data = [4, 3, 0, 5, 9, 0, 1, 0, 0, 0]
        self.assertEqual(perc_true(data), 0.5)

    def test_mean(self):
        data = [4, 3, 0, 5, 9, 0, 1, 0, 0, 0]
        self.assertEqual(mean(data), 2.2)

    def test_median(self):
        data = [4, 3, 0, 5, 9, 0, 1, 0, 0, 0]
        self.assertEqual(median(data), 0.5)

if __name__ == "__main__":
    unittest.main()
