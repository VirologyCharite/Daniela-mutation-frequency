#!/usr/bin/env python3

from dark.fasta import FastaReads
from math import log10
from heterogeneous_sequences import heterogeneousSites, compareToRef
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--homogenCutoff", type=float, default=1.0,
                    help="Cutoff fraction defining a site as homogeneous. 0.9 means 90\% or more is considered to be homogeneous.")
parser.add_argument("--reference",
                    help="Filename of the reference fasta sequence.")
args = parser.parse_args()

if 0 < args.homogenCutoff < 1:
    raise ValueError('--homogenCutoff needs to be between 0 and 1.')

reads = list(FastaReads('Alignment-CavV-strains.fasta'))

if len(set(len(read) for read in reads)) != 1:
    raise ValueError('Not all read lengths are the same.....')
else:
    length = len(reads[0]) 

if args.reference is not None:
    # A reference read has been provided.
    referenceread = list(FastaReads(args.reference))[0]
    mutationsdict = compareToRef(referenceread, reads, length)
    # C{defaultdict} containing mutation counts.
    totalmutations = sum(mutationsdict.values())
    mutationsSorted = sorted(zip(mutationsdict.values(), mutationsdict.keys()), reverse=True)
    # Mutation counts sorted by count.

print('Total count of mutations in sequences as compared to the reference sequence:')
for mutation in mutationsSorted:
    countofmutation = mutation[0]
    frequency = countofmutation/totalmutations
    print('%s mutations occur with a frequency of %f, %*d of %d.' % (mutation[1], frequency, 4, countofmutation, totalmutations))     

# Call heterogeneousSites function.
heterogeneous = heterogeneousSites(reads, length, args.homogenCutoff)

print('List of heterogeneous sites:')
width = int(log10(length)) + 1
for index, counts in heterogeneous.items():
    print('Site %*d has nucleotides %s' % (width, index + 1, str(counts)[9:-2]))

print('The number of heterogeneous sites is %d.' % len(heterogeneous))
