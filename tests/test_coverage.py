"""
This module contains tests for coverage methods.
"""

import unittest

from bp_quant.coverage import perc_true


class TestCoverage(unittest.TestCase):
    """Provides unit tests for coverage module."""

    def test_perc_true(self):
        """Test perc_true."""
        data = [4, 3, 0, 5, 9, 0, 1, 0, 0, 0]
        self.assertEqual(perc_true(data), 0.5)
