#!/usr/bin/env python

from argparse import ArgumentParser
import gzip
from random import choice, randint
import sys

def generate_random_sequence(length):
    DNA=""
    for count in range(length):
        DNA+=choice("CGTA")
    return DNA


def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))

class RandomSequenceGenerator(object):
    def __init__(self, seq_len, intervals, seq_count, read_len, read_count, min_insert_size):
        self.seq_len = seq_len
        self.intervals = intervals
        self.seq_count = seq_count
        self.read_len = read_len
        self.read_count = read_count
        self.min_insert_size = min_insert_size

    def run(self, fq1, fq2, tsv, no_gzip):
        r1_fq = None
        if no_gzip:
            r1_fq = open(fq1, "w")
        else:
            r1_fq = gzip.open(fq1, "wt")

        r2_fq = None
        if fq2:
            if no_gzip:
                r2_fq = open(fq2, "w")
            else:
                r2_fq = gzip.open(fq2, "wt")

        with open(tsv, "w") as outf:
            outf.write("name\tsequence\tposition\n")
            for i in range(self.seq_count):
                dna = generate_random_sequence(self.seq_len)


                outf.write("seq{}\t{}\t{}\n".format(i+1, dna, self.intervals))


                max_range = len(dna) - self.read_len
                if fq2:
                    max_range = len(dna) - self.read_len*2 - self.min_insert_size

                n = 0
                while n < self.read_count:
                    r1_start = randint(0, max_range)
                    r1_stop = r1_start + self.read_len
                    
                    r2_start = randint(r1_stop + self.min_insert_size, len(dna) - self.read_len)
                    r2_stop = r2_start + self.read_len

                    r1_seq = dna[r1_start:r1_stop]
                    r1_qual = "".join([choice("ABCDEFGHI") for x in range(self.read_len)])
                    r2_seq = rev_compl(dna[r2_start:r2_stop])
                    r2_qual = "".join([choice("ABCDEFGHI") for x in range(self.read_len)])
                    #r2_seq = dna[r2_start:r2_stop]
                    n += 1

                    #@HISEQ:218:C31BAACXX:1:1101:1380:2232 2:N:0:TGACCA
                    #GCGGGCAATGAGGAGAACTACGGTGTGAAGACTACCTATGATAGCAGTCT
                    #+
                    #@CCFFFFFFHFFHFD?FFHGFGBBBHEGGI@FHGIIIGIG@HEHIGHDGH

                    r1_fq.write("@HISEQ:{} 1:N:0:ACGTAC\n{}\n+\n{}\n".format(n, r1_seq, r1_qual))
        
                    if r2_fq:
                        r2_fq.write("@HISEQ:{} 2:N:0:ACGTAC\n{}\n+\n{}\n".format(n, r2_seq, r2_qual))

        r1_fq.close()

        if r2_fq:
            r2_fq.close()    

def main():

    parser = ArgumentParser(description="Generates random reads based on random sequence")
    parser.add_argument("-c", "--seq_len", dest="seq_len", type=int, help="Specify sequence length", default=1000)
    parser.add_argument("-p", "--intervals", dest="intervals", help="Specify comma-separated list of points of interest", default="0,200,400")
    parser.add_argument("-s", "--seq_count", dest="seq_count", type=int, help="Specify number of sequences to generate", default=10)
    parser.add_argument("-i", "--min_insert_size", dest="min_insert_size", type=int, help="Specify insert size for paired end FASTQs", default=0)
    parser.add_argument("-r", "--read_len", dest="read_len", type=int, help="Specify read length of FASTQ reads", default=50)
    parser.add_argument("-n", "--read_count", dest="read_count", type=int, help="Specify number of reads to generate", default=1000)
    parser.add_argument("-t", "--tsv", dest="tsv", help="Specify TSV output table path", required=True)
    parser.add_argument("-1", "--fastq1", dest="fq1", help="Specify path to FQ1", required=True)
    parser.add_argument("-2", "--fastq2", dest="fq2", help="Specify path to FQ2")
    parser.add_argument("--no_gzip", dest="no_gzip", action="store_true", help="Specify if output shall not be gzip compressed")

    args = parser.parse_args()

    rsg = RandomSequenceGenerator(args.seq_len, args.intervals, args.seq_count, args.read_len, args.read_count, args.min_insert_size)

    rsg.run(args.fq1, args.fq2, args.tsv, args.no_gzip)



if __name__ == "__main__":
    main()
