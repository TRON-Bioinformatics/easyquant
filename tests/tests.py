#!/usr/bin/env python

import os
import unittest

SEQ_TABLE_FILE = os.path.join("example_data", "CLDN18_Context_seq.csv")

from easy_quant.requantify import perc_true, mean, median, get_seq_to_pos, classify_read


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

        
    def test_get_seq_to_pos(self):
        result = {
            'CLDN18_1': [
                ('0_400', 0, 400),
                ('400_786', 400, 786)
            ],
            'CLDN18_2': [
                ('0_361', 0, 361),
                ('361_747', 361, 747)
            ],
            'CLDN18_total': [
                ('0_400', 0, 400),
                ('400_786', 400, 786)
            ],
            'CLDN18_1_fake': [
                ('0_400', 0, 400),
                ('400_786', 400, 786)
            ],
            'CLDN18_2_fake': [
                ('0_361', 0, 361),
                ('361_747', 361, 747)
            ],
            'HPRT1':[
                ('0_400', 0, 400),
                ('400_793', 400, 793)
            ]
        }
        self.assertEqual(get_seq_to_pos(SEQ_TABLE_FILE), result)

        
    def test_classify_read(self):
        aln_pairs = [(0, 365, 'C'), (1, 366, 'C'), (2, 367, 'C'), (3, 368, 'T'), (4, 369, 'A'), (5, 370, 'T'), (6, 371, 'T'), (7, 372, 'T'), (8, 373, 'C'), (9, 374, 'A'), (10, 375, 'C'), (11, 376, 'C'), (12, 377, 'A'), (13, 378, 'T'), (14, 379, 'C'), (15, 380, 'C'), (16, 381, 'T'), (17, 382, 'G'), (18, 383, 'G'), (19, 384, 'G'), (20, 385, 'A'), (21, 386, 'C'), (22, 387, 'T'), (23, 388, 'T'), (24, 389, 'C'), (25, 390, 'C'), (26, 391, 'A'), (27, 392, 'G'), (28, 393, 'C'), (29, 394, 'C'), (30, 395, 'A'), (31, 396, 'T'), (32, 397, 'G'), (33, 398, 'C'), (34, 399, 'T'), (35, 400, 'G'), (36, 401, 'C'), (37, 402, 'A'), (38, 403, 'G'), (39, 404, 'G'), (40, 405, 'C'), (41, 406, 'A'), (42, 407, 'G'), (43, 408, 'T'), (44, 409, 'G'), (45, 410, 'C'), (46, 411, 'G'), (47, 412, 'A'), (48, 413, 'G'), (49, 414, 'C'), (50, 415, 'C')]
        interval = [
            ('0_400', 0, 400),
            ('400_786', 400, 786)
        ]
        result = {"junc": True, "within": False, "interval": "0_400", "anchor": 15}
        self.assertEqual(classify_read(365, 415, aln_pairs, interval, True, 10), result)


        aln_pairs = [(0, 609, 'C'), (1, 610, 'T'), (2, 611, 'C'), (3, 612, 'A'), (4, 613, 'C'), (5, 614, 'C'), (6, 615, 'C'), (7, 616, 'A'), (8, 617, 'A'), (9, 618, 'A'), (10, 619, 'A'), (11, 620, 'A'), (12, 621, 'A'), (13, None, None), (14, 622, 'C'), (15, 623, 'A'), (16, 624, 'A'), (17, 625, 'G'), (18, 626, 'G'), (19, 627, 'A'), (20, 628, 'G'), (21, 629, 'A'), (22, 630, 'T'), (23, 631, 'C'), (24, 632, 'C'), (25, 633, 'C'), (26, 634, 'A'), (27, 635, 'T'), (28, 636, 'C'), (29, 637, 'T'), (30, 638, 'A'), (31, 639, 'G'), (32, 640, 'A'), (33, 641, 'T'), (34, 642, 'T'), (35, 643, 'T'), (36, 644, 'C'), (37, 645, 'T'), (38, 646, 'T'), (39, 647, 'C'), (40, 648, 'T'), (41, 649, 'T'), (42, 650, 'G'), (43, 651, 'C'), (44, 652, 'T'), (45, 653, 'T'), (46, 654, 'T'), (47, 655, 'T'), (48, 656, 'G'), (49, 657, 'A'), (50, 658, 'C')]

        interval = [
            ('0_400', 0, 400),
            ('400_786', 400, 786)
        ]
        result = {"junc": False, "within": True, "interval": "400_786", "anchor": 0}
        self.assertEqual(classify_read(609, 658, aln_pairs, interval, True, 10), result)


        # GTTGGACTATGCATGATAAGCAAGGTGAAGTAAGACTCAAATGTCTTACT
        # ATCATTTCGACATACAAGCACCCTGGCAGCTTCAGGAAAACCAGGTGGAA


        # TGGAGATTATCCACTTACCATGGCTGGTCCTCAGTGGAAGAAGTTCAAAT
        # CCAGTTTTTGTGAATTCATTGGCGTGTTAGTACGGCAATGTCAATATAGT
        # ATCATATATGATGAGTATATGATGGATACAGTCATTTCACTTCTTACAGG
        # ATTGTCTGACTCACAAGTCAG    AGCATTTCGACATACAAGCACCCTGGCAGCTTCAGGAAAATCAAGATGAAA
        # TAGAAAATATGATGAATGCAATATTTAAAGGAGTGTTTGTACATAGATACCGTGATGCGATAGCTGAAATTCGAGCTATTTGCATTGAAGAGATTGGCATTTGGATGAAGATGTATAGTGATGCCTTTCTTAATGACAGTTATTTAAAATATGTTGGTTGGACTATGCATGATAAGCA

        # 172
        aln_pairs = [
            (0, 172, 'A'),
            (1, 173, 'T'),
            (2, 174, 'C'),
            (3, 175, 'A'),
            (4, 176, 'T'),
            (5, 177, 'T'),
            (6, 178, 'T'),
            (7, 179, 'C'),
            (8, 180, 'G'),
            (9, 181, 'A'),
            (10, 182, 'C'),
            (11, 183, 'A'),
            (12, 184, 'T'),
            (13, 185, 'A'),
            (14, 186, 'C'),
            (15, 187, 'A'),
            (16, 188, 'A'),
            (17, 189, 'G'),
            (18, 190, 'C'),
            (19, 191, 'A'),
            (20, 192, 'C'),
            (21, 193, 'C'),
            (22, 194, 'C'),
            (23, 195, 'T'),
            (24, 196, 'G'),
            (25, 197, 'G'),
            (26, 198, 'C'),
            (27, 199, 'A'),
            (28, 200, 'G'),
            (29, 201, 'C'),
            (30, 202, 'T'),
            (31, 203, 'T'),
            (32, 204, 'C'),
            (33, 205, 'A'),
            (34, 206, 'G'),
            (35, 207, 'G'),
            (36, 208, 'A'),
            (37, 209, 'A'),
            (38, 210, 'A'),
            (39, 211, 'T'),
            (40, 212, 'C'),
            (41, 213, 'A'),
            (42, 214, 'A'),
            (43, 215, 'G'),
            (44, 216, 'A'),
            (45, 217, 'T'),  # TCAAGATGAAA
            (46, 218, 'G'),
            (47, 219, 'A'),
            (48, 220, 'A'),
            (49, 221, 'A')
        ]

        interval = [
            ('0_200', 0, 200),
            ('200_400', 200, 400)
        ]
        result = {"junc": True, "interval": "0_200", "within": False, "anchor": 22}
        self.assertEqual(classify_read(172, 222, aln_pairs, interval, True, 10), result)
        

if __name__ == "__main__":
    unittest.main()
