#!/usr/bin/env python

import os
import unittest

SEQ_TABLE_FILE = os.path.join("example_data", "CLDN18_Context_seq.csv")
BAM_FILE = os.path.join("example_data", "example_rna-seq.bam")

from bp_quant.requantify import Quantification, process_secondary_alignments, get_aligner, perc_true, get_seq_to_pos, classify_read


class TestRequantify(unittest.TestCase):

    def test_get_aligner(self):
        self.assertEqual(get_aligner(BAM_FILE), "star")


    def test_perc_true(self):
        data = [4, 3, 0, 5, 9, 0, 1, 0, 0, 0]
        self.assertEqual(perc_true(data), 0.5)


    def test_process_secondary_alignments(self):
        read_dict = {
            '1ef0b6243c3dcef97f849fd46d20ac25': {
                'A00259:290:HKKM5DMXY:1:1240:20383:22639': {
                    'R1': [{
                        'query_name': 'A00259:290:HKKM5DMXY:1:1240:20383:22639', 
                        'first_in_pair': True, 
                        'unmapped': False, 
                        'reference_name': '1ef0b6243c3dcef97f849fd46d20ac25', 
                        'flag': 355, 
                        'cigar': '100M', 
                        'start': 141, 
                        'stop': 241, 
                        'pairs': [(0, 141, 'C'), (1, 142, 'G'), (2, 143, 'G'), (3, 144, 'G'), (4, 145, 'G'), (5, 146, 'A'), (6, 147, 'G'), (7, 148, 'G'), (8, 149, 'T'), (9, 150, 'C'), (10, 151, 'G'), (11, 152, 'C'), (12, 153, 'C'), (13, 154, 'C'), (14, 155, 'C'), (15, 156, 'C'), (16, 157, 'G'), (17, 158, 'G'), (18, 159, 'G'), (19, 160, 'C'), (20, 161, 'G'), (21, 162, 'G'), (22, 163, 'C'), (23, 164, 'G'), (24, 165, 'A'), (25, 166, 'A'), (26, 167, 'G'), (27, 168, 'C'), (28, 169, 'a'), (29, 170, 'G'), (30, 171, 'C'), (31, 172, 'C'), (32, 173, 'A'), (33, 174, 'T'), (34, 175, 'T'), (35, 176, 'A'), (36, 177, 'T'), (37, 178, 'G'), (38, 179, 'C'), (39, 180, 'A'), (40, 181, 'G'), (41, 182, 'G'), (42, 183, 'C'), (43, 184, 'G'), (44, 185, 'A'), (45, 186, 'A'), (46, 187, 'T'), (47, 188, 'G'), (48, 189, 'A'), (49, 190, 'G'), (50, 191, 'G'), (51, 192, 'G'), (52, 193, 'A'), (53, 194, 'A'), (54, 195, 'T'), (55, 196, 'T'), (56, 197, 'G'), (57, 198, 'A'), (58, 199, 'G'), (59, 200, 'G'), (60, 201, 'T'), (61, 202, 'T'), (62, 203, 'C'), (63, 204, 'G'), (64, 205, 'T'), (65, 206, 'C'), (66, 207, 'T'), (67, 208, 'A'), (68, 209, 'G'), (69, 210, 'T'), (70, 211, 'A'), (71, 212, 'A'), (72, 213, 'A'), (73, 214, 'G'), (74, 215, 'A'), (75, 216, 'G'), (76, 217, 'G'), (77, 218, 'A'), (78, 219, 'G'), (79, 220, 'A'), (80, 221, 'G'), (81, 222, 'A'), (82, 223, 'A'), (83, 224, 'T'), (84, 225, 'A'), (85, 226, 'T'), (86, 227, 'G'), (87, 228, 'A'), (88, 229, 'T'), (89, 230, 'G'), (90, 231, 'A'), (91, 232, 'A'), (92, 233, 'G'), (93, 234, 'A'), (94, 235, 'G'), (95, 236, 'A'), (96, 237, 'C'), (97, 238, 'T'), (98, 239, 'G'), (99, 240, 'T')]
                    }], 
                    'R2': [{
                        'query_name': 'A00259:290:HKKM5DMXY:1:1240:20383:22639', 
                        'first_in_pair': False, 
                        'unmapped': False, 
                        'reference_name': '1ef0b6243c3dcef97f849fd46d20ac25', 
                        'flag': 403, 
                        'cigar': '100M', 
                        'start': 233, 
                        'stop': 333, 
                        'pairs': [(0, 233, 'G'), (1, 234, 'A'), (2, 235, 'G'), (3, 236, 'A'), (4, 237, 'C'), (5, 238, 'T'), (6, 239, 'G'), (7, 240, 'T'), (8, 241, 'T'), (9, 242, 'C'), (10, 243, 'G'), (11, 244, 'A'), (12, 245, 'T'), (13, 246, 'G'), (14, 247, 'G'), (15, 248, 'G'), (16, 249, 'C'), (17, 250, 'A'), (18, 251, 'G'), (19, 252, 'A'), (20, 253, 'T'), (21, 254, 'G'), (22, 255, 'C'), (23, 256, 'T'), (24, 257, 'G'), (25, 258, 'T'), (26, 259, 'C'), (27, 260, 'A'), (28, 261, 'T'), (29, 262, 'A'), (30, 263, 'G'), (31, 264, 'C'), (32, 265, 'T'), (33, 266, 'G'), (34, 267, 'C'), (35, 268, 'A'), (36, 269, 'G'), (37, 270, 'G'), (38, 271, 'A'), (39, 272, 'G'), (40, 273, 'G'), (41, 274, 'T'), (42, 275, 'G'), (43, 276, 'A'), (44, 277, 'T'), (45, 278, 'G'), (46, 279, 'G'), (47, 280, 'C'), (48, 281, 'A'), (49, 282, 'C'), (50, 283, 'A'), (51, 284, 'A'), (52, 285, 'T'), (53, 286, 'G'), (54, 287, 'C'), (55, 288, 'T'), (56, 289, 'G'), (57, 290, 'C'), (58, 291, 'T'), (59, 292, 'G'), (60, 293, 'G'), (61, 294, 'C'), (62, 295, 'A'), (63, 296, 'G'), (64, 297, 'C'), (65, 298, 'G'), (66, 299, 'A'), (67, 300, 'G'), (68, 301, 'T'), (69, 302, 'A'), (70, 303, 'A'), (71, 304, 'A'), (72, 305, 'G'), (73, 306, 'T'), (74, 307, 'C'), (75, 308, 'T'), (76, 309, 'T'), (77, 310, 'G'), (78, 311, 'G'), (79, 312, 'A'), (80, 313, 'C'), (81, 314, 'A'), (82, 315, 'G'), (83, 316, 'A'), (84, 317, 'C'), (85, 318, 'T'), (86, 319, 'T'), (87, 320, 'A'), (88, 321, 'A'), (89, 322, 'A'), (90, 323, 'C'), (91, 324, 'C'), (92, 325, 'A'), (93, 326, 'G'), (94, 327, 'T'), (95, 328, 'T'), (96, 329, 'A'), (97, 330, 'T'), (98, 331, 'A'), (99, 332, 'G')]
                    }]
                }
            }, 
            '49835f1bb383b6c7f76a2099892fa4ff': {
                'A00259:290:HKKM5DMXY:2:2332:4164:35070': {
                    'R1': [], 
                    'R2': [{
                        'query_name': 'A00259:290:HKKM5DMXY:2:2332:4164:35070', 
                        'first_in_pair': False, 
                        'unmapped': False, 
                        'reference_name': '49835f1bb383b6c7f76a2099892fa4ff', 
                        'flag': 401, 
                        'cigar': '100M', 
                        'start': 284, 
                        'stop': 384, 
                        'pairs': [(0, 284, 'C'), (1, 285, 'C'), (2, 286, 'T'), (3, 287, 'G'), (4, 288, 'C'), (5, 289, 'G'), (6, 290, 'T'), (7, 291, 'C'), (8, 292, 'T'), (9, 293, 'G'), (10, 294, 'T'), (11, 295, 'G'), (12, 296, 'C'), (13, 297, 'T'), (14, 298, 'C'), (15, 299, 'C'), (16, 300, 'G'), (17, 301, 'A'), (18, 302, 'G'), (19, 303, 'G'), (20, 304, 'T'), (21, 305, 'G'), (22, 306, 'C'), (23, 307, 'C'), (24, 308, 'C'), (25, 309, 'A'), (26, 310, 'T'), (27, 311, 'T'), (28, 312, 't'), (29, 313, 'T'), (30, 314, 'G'), (31, 315, 'G'), (32, 316, 'A'), (33, 317, 'G'), (34, 318, 'G'), (35, 319, 'A'), (36, 320, 'C'), (37, 321, 'A'), (38, 322, 'C'), (39, 323, 'G'), (40, 324, 'C'), (41, 325, 'T'), (42, 326, 'G'), (43, 327, 'A'), (44, 328, 'T'), (45, 329, 'G'), (46, 330, 'C'), (47, 331, 'G'), (48, 332, 'C'), (49, 333, 'A'), (50, 334, 'T'), (51, 335, 'C'), (52, 336, 'C'), (53, 337, 'T'), (54, 338, 'G'), (55, 339, 'G'), (56, 340, 'T'), (57, 341, 'C'), (58, 342, 'A'), (59, 343, 'T'), (60, 344, 'C'), (61, 345, 'G'), (62, 346, 'G'), (63, 347, 'G'), (64, 348, 'C'), (65, 349, 'T'), (66, 350, 'G'), (67, 351, 'T'), (68, 352, 'C'), (69, 353, 'C'), (70, 354, 'C'), (71, 355, 'G'), (72, 356, 'G'), (73, 357, 'G'), (74, 358, 'A'), (75, 359, 'G'), (76, 360, 'C'), (77, 361, 'T'), (78, 362, 'C'), (79, 363, 'C'), (80, 364, 'C'), (81, 365, 'G'), (82, 366, 'C'), (83, 367, 'T'), (84, 368, 'C'), (85, 369, 'G'), (86, 370, 'G'), (87, 371, 'G'), (88, 372, 'C'), (89, 373, 'C'), (90, 374, 'T'), (91, 375, 'G'), (92, 376, 'C'), (93, 377, 'G'), (94, 378, 'G'), (95, 379, 'A'), (96, 380, 'C'), (97, 381, 'G'), (98, 382, 'C'), (99, 383, 'C')]
                    }]
                }
            }, 
            '38f796f6fe0f1c7b891c94043ef94e33': {
                'A00259:290:HKKM5DMXY:2:2332:4164:35070': {
                    'R1': [{
                        'query_name': 'A00259:290:HKKM5DMXY:2:2332:4164:35070', 
                        'first_in_pair': True, 
                        'unmapped': False, 
                        'reference_name': '38f796f6fe0f1c7b891c94043ef94e33', 
                        'flag': 353, 
                        'cigar': '100M', 
                        'start': 23, 
                        'stop': 123, 
                        'pairs': [(0, 23, 'G'), (1, 24, 'C'), (2, 25, 'C'), (3, 26, 'G'), (4, 27, 'G'), (5, 28, 'C'), (6, 29, 'A'), (7, 30, 'T'), (8, 31, 'C'), (9, 32, 'G'), (10, 33, 'C'), (11, 34, 'C'), (12, 35, 'T'), (13, 36, 'G'), (14, 37, 'G'), (15, 38, 'G'), (16, 39, 'A'), (17, 40, 'C'), (18, 41, 'A'), (19, 42, 'A'), (20, 43, 'A'), (21, 44, 'G'), (22, 45, 'G'), (23, 46, 'A'), (24, 47, 'G'), (25, 48, 'A'), (26, 49, 'A'), (27, 50, 'A'), (28, 51, 'A'), (29, 52, 'G'), (30, 53, 'A'), (31, 54, 'G'), (32, 55, 'G'), (33, 56, 'A'), (34, 57, 'A'), (35, 58, 'A'), (36, 59, 'C'), (37, 60, 'C'), (38, 61, 'A'), (39, 62, 'G'), (40, 63, 'A'), (41, 64, 'T'), (42, 65, 'T'), (43, 66, 'G'), (44, 67, 'C'), (45, 68, 'C'), (46, 69, 'G'), (47, 70, 'C'), (48, 71, 'C'), (49, 72, 'A'), (50, 73, 'T'), (51, 74, 'C'), (52, 75, 'C'), (53, 76, 'A'), (54, 77, 'G'), (55, 78, 'C'), (56, 79, 'G'), (57, 80, 'G'), (58, 81, 'G'), (59, 82, 'A'), (60, 83, 'T'), (61, 84, 'G'), (62, 85, 'C'), (63, 86, 'C'), (64, 87, 'G'), (65, 88, 'T'), (66, 89, 'C'), (67, 90, 'T'), (68, 91, 'G'), (69, 92, 'G'), (70, 93, 'T'), (71, 94, 'G'), (72, 95, 'G'), (73, 96, 'C'), (74, 97, 'T'), (75, 98, 'C'), (76, 99, 'C'), (77, 100, 'A'), (78, 101, 'C'), (79, 102, 'A'), (80, 103, 'C'), (81, 104, 'T'), (82, 105, 'G'), (83, 106, 'T'), (84, 107, 'G'), (85, 108, 'G'), (86, 109, 'T'), (87, 110, 'C'), (88, 111, 'C'), (89, 112, 'C'), (90, 113, 'C'), (91, 114, 'T'), (92, 115, 'C'), (93, 116, 'C'), (94, 117, 'A'), (95, 118, 'T'), (96, 119, 'C'), (97, 120, 'A'), (98, 121, 'G'), (99, 122, 'C')]
                    }], 
                    'R2': []
                }
            }
        }
        read_pairings = [
            ({
                'query_name': 'A00259:290:HKKM5DMXY:1:1240:20383:22639', 
                'first_in_pair': True, 
                'unmapped': False, 
                'reference_name': '1ef0b6243c3dcef97f849fd46d20ac25', 
                'flag': 355, 
                'cigar': '100M', 
                'start': 141, 
                'stop': 241, 
                'pairs': [(0, 141, 'C'), (1, 142, 'G'), (2, 143, 'G'), (3, 144, 'G'), (4, 145, 'G'), (5, 146, 'A'), (6, 147, 'G'), (7, 148, 'G'), (8, 149, 'T'), (9, 150, 'C'), (10, 151, 'G'), (11, 152, 'C'), (12, 153, 'C'), (13, 154, 'C'), (14, 155, 'C'), (15, 156, 'C'), (16, 157, 'G'), (17, 158, 'G'), (18, 159, 'G'), (19, 160, 'C'), (20, 161, 'G'), (21, 162, 'G'), (22, 163, 'C'), (23, 164, 'G'), (24, 165, 'A'), (25, 166, 'A'), (26, 167, 'G'), (27, 168, 'C'), (28, 169, 'a'), (29, 170, 'G'), (30, 171, 'C'), (31, 172, 'C'), (32, 173, 'A'), (33, 174, 'T'), (34, 175, 'T'), (35, 176, 'A'), (36, 177, 'T'), (37, 178, 'G'), (38, 179, 'C'), (39, 180, 'A'), (40, 181, 'G'), (41, 182, 'G'), (42, 183, 'C'), (43, 184, 'G'), (44, 185, 'A'), (45, 186, 'A'), (46, 187, 'T'), (47, 188, 'G'), (48, 189, 'A'), (49, 190, 'G'), (50, 191, 'G'), (51, 192, 'G'), (52, 193, 'A'), (53, 194, 'A'), (54, 195, 'T'), (55, 196, 'T'), (56, 197, 'G'), (57, 198, 'A'), (58, 199, 'G'), (59, 200, 'G'), (60, 201, 'T'), (61, 202, 'T'), (62, 203, 'C'), (63, 204, 'G'), (64, 205, 'T'), (65, 206, 'C'), (66, 207, 'T'), (67, 208, 'A'), (68, 209, 'G'), (69, 210, 'T'), (70, 211, 'A'), (71, 212, 'A'), (72, 213, 'A'), (73, 214, 'G'), (74, 215, 'A'), (75, 216, 'G'), (76, 217, 'G'), (77, 218, 'A'), (78, 219, 'G'), (79, 220, 'A'), (80, 221, 'G'), (81, 222, 'A'), (82, 223, 'A'), (83, 224, 'T'), (84, 225, 'A'), (85, 226, 'T'), (86, 227, 'G'), (87, 228, 'A'), (88, 229, 'T'), (89, 230, 'G'), (90, 231, 'A'), (91, 232, 'A'), (92, 233, 'G'), (93, 234, 'A'), (94, 235, 'G'), (95, 236, 'A'), (96, 237, 'C'), (97, 238, 'T'), (98, 239, 'G'), (99, 240, 'T')]
            },
            {
                'query_name': 'A00259:290:HKKM5DMXY:1:1240:20383:22639',
                'first_in_pair': False,
                'unmapped': False,
                'reference_name': '1ef0b6243c3dcef97f849fd46d20ac25',
                'flag': 403,
                'cigar': '100M',
                'start': 233,
                'stop': 333,
                'pairs': [(0, 233, 'G'), (1, 234, 'A'), (2, 235, 'G'), (3, 236, 'A'), (4, 237, 'C'), (5, 238, 'T'), (6, 239, 'G'), (7, 240, 'T'), (8, 241, 'T'), (9, 242, 'C'), (10, 243, 'G'), (11, 244, 'A'), (12, 245, 'T'), (13, 246, 'G'), (14, 247, 'G'), (15, 248, 'G'), (16, 249, 'C'), (17, 250, 'A'), (18, 251, 'G'), (19, 252, 'A'), (20, 253, 'T'), (21, 254, 'G'), (22, 255, 'C'), (23, 256, 'T'), (24, 257, 'G'), (25, 258, 'T'), (26, 259, 'C'), (27, 260, 'A'), (28, 261, 'T'), (29, 262, 'A'), (30, 263, 'G'), (31, 264, 'C'), (32, 265, 'T'), (33, 266, 'G'), (34, 267, 'C'), (35, 268, 'A'), (36, 269, 'G'), (37, 270, 'G'), (38, 271, 'A'), (39, 272, 'G'), (40, 273, 'G'), (41, 274, 'T'), (42, 275, 'G'), (43, 276, 'A'), (44, 277, 'T'), (45, 278, 'G'), (46, 279, 'G'), (47, 280, 'C'), (48, 281, 'A'), (49, 282, 'C'), (50, 283, 'A'), (51, 284, 'A'), (52, 285, 'T'), (53, 286, 'G'), (54, 287, 'C'), (55, 288, 'T'), (56, 289, 'G'), (57, 290, 'C'), (58, 291, 'T'), (59, 292, 'G'), (60, 293, 'G'), (61, 294, 'C'), (62, 295, 'A'), (63, 296, 'G'), (64, 297, 'C'), (65, 298, 'G'), (66, 299, 'A'), (67, 300, 'G'), (68, 301, 'T'), (69, 302, 'A'), (70, 303, 'A'), (71, 304, 'A'), (72, 305, 'G'), (73, 306, 'T'), (74, 307, 'C'), (75, 308, 'T'), (76, 309, 'T'), (77, 310, 'G'), (78, 311, 'G'), (79, 312, 'A'), (80, 313, 'C'), (81, 314, 'A'), (82, 315, 'G'), (83, 316, 'A'), (84, 317, 'C'), (85, 318, 'T'), (86, 319, 'T'), (87, 320, 'A'), (88, 321, 'A'), (89, 322, 'A'), (90, 323, 'C'), (91, 324, 'C'), (92, 325, 'A'), (93, 326, 'G'), (94, 327, 'T'), (95, 328, 'T'), (96, 329, 'A'), (97, 330, 'T'), (98, 331, 'A'), (99, 332, 'G')]
            }),
            ({
                "reference_name": "49835f1bb383b6c7f76a2099892fa4ff",
                "query_name": "A00259:290:HKKM5DMXY:2:2332:4164:35070",
                "unmapped": True,
                "flag": 325,
                "start": -1,
                "stop": -1,
                "pairs": None,
                "cigar": None
            },
            {
                'query_name': 'A00259:290:HKKM5DMXY:2:2332:4164:35070', 
                'first_in_pair': False, 
                'unmapped': False, 
                'reference_name': '49835f1bb383b6c7f76a2099892fa4ff', 
                'flag': 401, 
                'cigar': '100M', 
                'start': 284, 
                'stop': 384, 
                'pairs': [(0, 284, 'C'), (1, 285, 'C'), (2, 286, 'T'), (3, 287, 'G'), (4, 288, 'C'), (5, 289, 'G'), (6, 290, 'T'), (7, 291, 'C'), (8, 292, 'T'), (9, 293, 'G'), (10, 294, 'T'), (11, 295, 'G'), (12, 296, 'C'), (13, 297, 'T'), (14, 298, 'C'), (15, 299, 'C'), (16, 300, 'G'), (17, 301, 'A'), (18, 302, 'G'), (19, 303, 'G'), (20, 304, 'T'), (21, 305, 'G'), (22, 306, 'C'), (23, 307, 'C'), (24, 308, 'C'), (25, 309, 'A'), (26, 310, 'T'), (27, 311, 'T'), (28, 312, 't'), (29, 313, 'T'), (30, 314, 'G'), (31, 315, 'G'), (32, 316, 'A'), (33, 317, 'G'), (34, 318, 'G'), (35, 319, 'A'), (36, 320, 'C'), (37, 321, 'A'), (38, 322, 'C'), (39, 323, 'G'), (40, 324, 'C'), (41, 325, 'T'), (42, 326, 'G'), (43, 327, 'A'), (44, 328, 'T'), (45, 329, 'G'), (46, 330, 'C'), (47, 331, 'G'), (48, 332, 'C'), (49, 333, 'A'), (50, 334, 'T'), (51, 335, 'C'), (52, 336, 'C'), (53, 337, 'T'), (54, 338, 'G'), (55, 339, 'G'), (56, 340, 'T'), (57, 341, 'C'), (58, 342, 'A'), (59, 343, 'T'), (60, 344, 'C'), (61, 345, 'G'), (62, 346, 'G'), (63, 347, 'G'), (64, 348, 'C'), (65, 349, 'T'), (66, 350, 'G'), (67, 351, 'T'), (68, 352, 'C'), (69, 353, 'C'), (70, 354, 'C'), (71, 355, 'G'), (72, 356, 'G'), (73, 357, 'G'), (74, 358, 'A'), (75, 359, 'G'), (76, 360, 'C'), (77, 361, 'T'), (78, 362, 'C'), (79, 363, 'C'), (80, 364, 'C'), (81, 365, 'G'), (82, 366, 'C'), (83, 367, 'T'), (84, 368, 'C'), (85, 369, 'G'), (86, 370, 'G'), (87, 371, 'G'), (88, 372, 'C'), (89, 373, 'C'), (90, 374, 'T'), (91, 375, 'G'), (92, 376, 'C'), (93, 377, 'G'), (94, 378, 'G'), (95, 379, 'A'), (96, 380, 'C'), (97, 381, 'G'), (98, 382, 'C'), (99, 383, 'C')]
            }),
            ({
                'query_name': 'A00259:290:HKKM5DMXY:2:2332:4164:35070', 
                'first_in_pair': True, 
                'unmapped': False, 
                'reference_name': '38f796f6fe0f1c7b891c94043ef94e33', 
                'flag': 353, 
                'cigar': '100M', 
                'start': 23, 
                'stop': 123, 
                'pairs': [(0, 23, 'G'), (1, 24, 'C'), (2, 25, 'C'), (3, 26, 'G'), (4, 27, 'G'), (5, 28, 'C'), (6, 29, 'A'), (7, 30, 'T'), (8, 31, 'C'), (9, 32, 'G'), (10, 33, 'C'), (11, 34, 'C'), (12, 35, 'T'), (13, 36, 'G'), (14, 37, 'G'), (15, 38, 'G'), (16, 39, 'A'), (17, 40, 'C'), (18, 41, 'A'), (19, 42, 'A'), (20, 43, 'A'), (21, 44, 'G'), (22, 45, 'G'), (23, 46, 'A'), (24, 47, 'G'), (25, 48, 'A'), (26, 49, 'A'), (27, 50, 'A'), (28, 51, 'A'), (29, 52, 'G'), (30, 53, 'A'), (31, 54, 'G'), (32, 55, 'G'), (33, 56, 'A'), (34, 57, 'A'), (35, 58, 'A'), (36, 59, 'C'), (37, 60, 'C'), (38, 61, 'A'), (39, 62, 'G'), (40, 63, 'A'), (41, 64, 'T'), (42, 65, 'T'), (43, 66, 'G'), (44, 67, 'C'), (45, 68, 'C'), (46, 69, 'G'), (47, 70, 'C'), (48, 71, 'C'), (49, 72, 'A'), (50, 73, 'T'), (51, 74, 'C'), (52, 75, 'C'), (53, 76, 'A'), (54, 77, 'G'), (55, 78, 'C'), (56, 79, 'G'), (57, 80, 'G'), (58, 81, 'G'), (59, 82, 'A'), (60, 83, 'T'), (61, 84, 'G'), (62, 85, 'C'), (63, 86, 'C'), (64, 87, 'G'), (65, 88, 'T'), (66, 89, 'C'), (67, 90, 'T'), (68, 91, 'G'), (69, 92, 'G'), (70, 93, 'T'), (71, 94, 'G'), (72, 95, 'G'), (73, 96, 'C'), (74, 97, 'T'), (75, 98, 'C'), (76, 99, 'C'), (77, 100, 'A'), (78, 101, 'C'), (79, 102, 'A'), (80, 103, 'C'), (81, 104, 'T'), (82, 105, 'G'), (83, 106, 'T'), (84, 107, 'G'), (85, 108, 'G'), (86, 109, 'T'), (87, 110, 'C'), (88, 111, 'C'), (89, 112, 'C'), (90, 113, 'C'), (91, 114, 'T'), (92, 115, 'C'), (93, 116, 'C'), (94, 117, 'A'), (95, 118, 'T'), (96, 119, 'C'), (97, 120, 'A'), (98, 121, 'G'), (99, 122, 'C')]
            },
            {
                "reference_name": "38f796f6fe0f1c7b891c94043ef94e33",
                "query_name": "A00259:290:HKKM5DMXY:2:2332:4164:35070",
                "unmapped": True,
                "flag": 389,
                "start": -1,
                "stop": -1,
                "pairs": None,
                "cigar": None
            })
        ]
        self.assertEqual(process_secondary_alignments(read_dict), read_pairings)
        

    def test_get_seq_to_pos(self):
        result = {
            'CLDN18_1': (
                [
                    ('0_400', 0, 400),
                    ('400_786', 400, 786)
                ], 
                "0,400,786"
            ),
            'CLDN18_2': (
                [
                    ('0_361', 0, 361),
                    ('361_747', 361, 747)
                ],
                "0,361,747"
            ),
            'CLDN18_total': (
                [
                    ('0_400', 0, 400),
                    ('400_786', 400, 786)
                ],
                "0,400,786"
            ),
            'CLDN18_1_fake': (
                [
                    ('0_400', 0, 400),
                    ('400_786', 400, 786)
                ],
                "0,400,786"
            ),
            'CLDN18_2_fake': (
                [
                    ('0_361', 0, 361),
                    ('361_747', 361, 747)
                ],
                "0,361,747"
            ),
            'HPRT1': (
                [
                    ('0_400', 0, 400),
                    ('400_793', 400, 793)
                ],
                "0,400,793"
            ),
            'HPRT1_dup': (
                [
                    ('0_400', 0, 400),
                    ('400_793', 400, 793)
                ],
                "0,400,793"
            ),
            'HPRT1_similar': (
                [
                    ('0_400', 0, 400),
                    ('400_793', 400, 793)
                ],
                "0,400,793"
            )
        }
        self.assertEqual(get_seq_to_pos(SEQ_TABLE_FILE), result)

        
    def test_classify_read(self):
        aln_pairs = [(0, 365, 'C'), (1, 366, 'C'), (2, 367, 'C'), (3, 368, 'T'), (4, 369, 'A'), (5, 370, 'T'), (6, 371, 'T'), (7, 372, 'T'), (8, 373, 'C'), (9, 374, 'A'), (10, 375, 'C'), (11, 376, 'C'), (12, 377, 'A'), (13, 378, 'T'), (14, 379, 'C'), (15, 380, 'C'), (16, 381, 'T'), (17, 382, 'G'), (18, 383, 'G'), (19, 384, 'G'), (20, 385, 'A'), (21, 386, 'C'), (22, 387, 'T'), (23, 388, 'T'), (24, 389, 'C'), (25, 390, 'C'), (26, 391, 'A'), (27, 392, 'G'), (28, 393, 'C'), (29, 394, 'C'), (30, 395, 'A'), (31, 396, 'T'), (32, 397, 'G'), (33, 398, 'C'), (34, 399, 'T'), (35, 400, 'G'), (36, 401, 'C'), (37, 402, 'A'), (38, 403, 'G'), (39, 404, 'G'), (40, 405, 'C'), (41, 406, 'A'), (42, 407, 'G'), (43, 408, 'T'), (44, 409, 'G'), (45, 410, 'C'), (46, 411, 'G'), (47, 412, 'A'), (48, 413, 'G'), (49, 414, 'C'), (50, 415, 'C')]
        interval = [
            ('0_400', 0, 400),
            ('400_786', 400, 786)
        ]
        result = {"junc": True, "within": False, "interval": "0_400", "anchor": 15, "nm": 0, "nm_in_bp_area": 0}
        self.assertEqual(classify_read(365, 415, aln_pairs, interval, True, 10), result)


        aln_pairs = [(0, 609, 'C'), (1, 610, 'T'), (2, 611, 'C'), (3, 612, 'A'), (4, 613, 'C'), (5, 614, 'C'), (6, 615, 'C'), (7, 616, 'A'), (8, 617, 'A'), (9, 618, 'A'), (10, 619, 'A'), (11, 620, 'A'), (12, 621, 'A'), (13, None, None), (14, 622, 'C'), (15, 623, 'A'), (16, 624, 'A'), (17, 625, 'G'), (18, 626, 'G'), (19, 627, 'A'), (20, 628, 'G'), (21, 629, 'A'), (22, 630, 'T'), (23, 631, 'C'), (24, 632, 'C'), (25, 633, 'C'), (26, 634, 'A'), (27, 635, 'T'), (28, 636, 'C'), (29, 637, 'T'), (30, 638, 'A'), (31, 639, 'G'), (32, 640, 'A'), (33, 641, 'T'), (34, 642, 'T'), (35, 643, 'T'), (36, 644, 'C'), (37, 645, 'T'), (38, 646, 'T'), (39, 647, 'C'), (40, 648, 'T'), (41, 649, 'T'), (42, 650, 'G'), (43, 651, 'C'), (44, 652, 'T'), (45, 653, 'T'), (46, 654, 'T'), (47, 655, 'T'), (48, 656, 'G'), (49, 657, 'A'), (50, 658, 'C')]

        interval = [
            ('0_400', 0, 400),
            ('400_786', 400, 786)
        ]
        result = {"junc": False, "within": True, "interval": "400_786", "anchor": 0, "nm": 0, "nm_in_bp_area": 0}
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
        result = {"junc": True, "interval": "0_200", "within": False, "anchor": 22, "nm": 0, "nm_in_bp_area": 0}
        self.assertEqual(classify_read(172, 222, aln_pairs, interval, True, 10), result)
        
    def test_classify_softjunc(self):
        aln_pairs = [(0, 609, 'C'), (1, 610, 'T'), (2, 611, 'C'), (3, 612, 'A'), (4, 613, 'C'), (5, 614, 'C'), (6, 615, 'C'), (7, 616, 'A'), (8, 617, 'A'), (9, 618, 'A'), (10, 619, 'A'), (11, 620, 'A'), (12, 621, 'A'), (13, None, None), (14, 622, 'C'), (15, 623, 'A'), (16, 624, 'A'), (17, 625, 'G'), (18, 626, 'G'), (19, 627, 'A'), (20, 628, 'G'), (21, 629, 'A'), (22, 630, 'T'), (23, 631, 'C'), (24, 632, 'C'), (25, 633, 'C'), (26, 634, 'A'), (27, 635, 'T'), (28, 636, 'C'), (29, 637, 'T'), (30, 638, 'A'), (31, 639, 'G'), (32, 640, 'A'), (33, 641, 'T'), (34, 642, 'T'), (35, 643, 'T'), (36, 644, 'C'), (37, 645, 'T'), (38, 646, 'T'), (39, 647, 'C'), (40, 648, 'T'), (41, 649, 'T'), (42, 650, 'G'), (43, 651, 'C'), (44, 652, 'T'), (45, 653, 'T'), (46, 654, 'T'), (47, 655, 'T'), (48, 656, 'G'), (49, 657, 'A'), (50, 658, 'C')]

        interval = [
            ('0_400', 0, 400),
            ('400_786', 400, 786)
        ]
        result = {"junc": False, "within": False, "interval": "", "anchor": 7, "nm": 0, "nm_in_bp_area": 0}
        self.assertEqual(classify_read(393, 444, aln_pairs, interval, True, 10), result)

        # J00128:23:H323KBBXX:6:2228:2735:12638   R1      339     HPRT1   0,400,793       393     444     51M     0       0
        # J00128:23:H323KBBXX:6:2228:2735:12638   R2      419     HPRT1   0,400,793       150     201     51M     0       0       within

if __name__ == "__main__":
    unittest.main()
