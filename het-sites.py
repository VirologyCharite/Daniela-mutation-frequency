#!/usr/bin/env python3

from __future__ import print_function, division

from collections import Counter
from dark.fasta import FastaReads
from math import log10
import argparse

reads = list(FastaReads('Alignment-CavV-strains.fasta'))

if len(set(len(read) for read in reads)) != 1:
    raise ValueError('Not all read lengths are the same.....')

length = len(reads[0])
width = int(log10(length)) + 1
total_heterogeneious = 0

for index in range(length):
    counts = Counter()
    for read in reads:
        counts[read.sequence[index]] += 1

    if len(counts) > 1:
        # The site is not homogenous = more than one nucleotide present.
        print('Site %*d has nucleotides %s' % (width, index + 1, counts))
        total_heterogeneious += 1

print('The total of heterogeneious sites is %i.' % (total_heterogeneious))

# add a parsing function, when to include a site as heterogeneious