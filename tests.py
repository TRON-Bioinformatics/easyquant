#!/usr/bin/env python

import unittest

from easy_quant import get_fastq_files, pair_fastq_files

class TestIOMethods(unittest.TestCase):

    def test_get_fastqs(self):
        pass

    def test_pair_fastqs(self):
        self.assertEqual(pair_fastq_files(["/a/b/c/test_R1.fastq.gz", "/a/b/c/test_R2.fastq.gz"]), (["/a/b/c/test_R1.fastq.gz"], ["/a/b/c/test_R2.fastq.gz"], ["test"]))

        self.assertEqual(pair_fastq_files(["/a/b/c/test.R1.fastp.fastq.gz", "/a/b/c/test.R2.fastp.fastq.gz"]), (["/a/b/c/test.R1.fastp.fastq.gz"], ["/a/b/c/test.R2.fastp.fastq.gz"], ["test"]))


if __name__ == "__main__":
    unittest.main()
