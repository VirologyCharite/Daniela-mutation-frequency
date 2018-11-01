#!/usr/bin/env python3

from __future__ import print_function, division

from collections import Counter
from dark.fasta import FastaReads
from math import log10
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--homogen", type=float, default=1,
                    help="cutoff fraction of base at site at which a site is considered to be homogenous.")
args = parser.parse_args()

if args.homogen > 1 or args.homogen < 0:
    raise ValueError('--homogen needs to be between 0 and 1.')

reads = list(FastaReads('Alignment-CavV-strains.fasta'))

if len(set(len(read) for read in reads)) != 1:
    raise ValueError('Not all read lengths are the same.....')

length = len(reads[0])
width = int(log10(length)) + 1
total_heterogeneious = 0

for index in range(length):
    counts = Counter()
    for read in reads:
        # Do not take gaps into account.
        if read.sequence[index] != '-':    
            counts[read.sequence[index]] += 1

    if len(counts) > 1:
        if (counts.most_common(1)[0][1])/sum(counts.values()) < args.homogen:
            # The site is not homogenous = more than one nucleotide present.       
            print('Site %*d has nucleotides %s' % (width, index + 1, str(counts)[9:-2]))
            total_heterogeneious += 1

print('The total of heterogeneious sites is %i.' % (total_heterogeneious))

