#!/usr/bin/env python3

from __future__ import print_function, division

from collections import Counter
from dark.fasta import FastaReads
from math import log10
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--homogen", type=float, default=1,
                    help="Cutoff fraction of base at site at which a site is considered to be homogenous.")
parser.add_argument("--ref", type=str, default="",
                    help="Filename of the reference fasta sequence.")
args = parser.parse_args()

if args.homogen > 1 or args.homogen < 0:
    raise ValueError('--homogen needs to be between 0 and 1.')

reads = list(FastaReads('Alignment-CavV-strains.fasta'))

if len(set(len(read) for read in reads)) != 1:
    raise ValueError('Not all read lengths are the same.....')

length = len(reads[0])
width = int(log10(length)) + 1

def hetsites():
    total_heterogeneious = 0

    for index in range(length):
        counts = Counter()

        for read in reads:
            if read.sequence[index] != '-':
            # Do not take gaps into account.   
                counts[read.sequence[index]] += 1

        if len(counts) > 1:
            # The site is not homogenous = more than one nucleotide present. 
            if (counts.most_common(1)[0][1])/sum(counts.values()) < args.homogen:
                # Most common base/total bases is the cutoff       
                print('Site %*d has nucleotides %s' % (width, index + 1, str(counts)[9:-2]))
                total_heterogeneious += 1

    print('The total of heterogeneious sites is %i.' % (total_heterogeneious))

def comparetoref():

    reference = list(FastaReads(args.ref))[0]

    #if len(reference.sequence) != length:
     #   raise ValueError('Reference read length is not the same as length of other reads.')

    totalmutations = 0
    totalsitesmutated = 0
    mutations = {
                'A to C' : 0,
                'A to G' : 0,
                'A to T' : 0,
                'C to A' : 0,
                'C to G' : 0,
                'C to T' : 0,
                'G to A' : 0,
                'G to C' : 0,
                'G to T' : 0,
                'T to A' : 0,
                'T to C' : 0,
                'T to G' : 0}

    for index in range(length):

        refbase = reference.sequence[index]
        # Exclude gaps.
        if refbase == '-':
            continue

        is_mutation = False

        for read in reads:
            base = read.sequence[index]

            if base == '-':
                # Exclude gaps.
                continue

            if base == refbase:
                continue

            else:
                is_mutation = True
                totalmutations += 1

            if refbase == 'A':
                if base == 'C':
                    mutations['A to C'] += 1
                elif base == 'G':
                    mutations['A to G'] += 1
                elif base == 'U' or base == 'T':
                    mutations['A to T'] += 1
            elif refbase == 'C':
                if base == 'A':
                    mutations['C to A'] += 1
                elif base == 'G':
                    mutations['C to G'] += 1
                elif base == 'U' or base == 'T':
                    mutations['C to T'] += 1
            elif refbase == 'G':
                if base == 'A':
                    mutations['G to A'] += 1
                elif base == 'C':
                    mutations['G to C'] += 1
                elif base == 'U' or base == 'T':
                    mutations['G to T'] += 1
            elif refbase == 'U' or refbase == 'T':
                if base == 'A':
                    mutations['T to A'] += 1
                elif base == 'C':
                    mutations['T to C'] += 1
                elif base == 'G':
                    mutations['T to G'] += 1
            
        if is_mutation == True:
            totalsitesmutated += 1

    allmutations = sum(mutations.values())
    if allmutations != totalmutations:
        # What kind of an error to raise here?
        raise ValueError('Allmutations and totalmutations should be the same.')

    print('At %f percent of sites, a mutation occurred.' % (totalsitesmutated/length))

    for mutation in mutations:
        numberofmut = mutations[mutation]
        frequency = numberofmut/totalmutations
        print('%s mutations occur with a frequency of %f, %*i of %i.' % (mutation, frequency, 4, numberofmut, totalmutations)) 

if args.ref != "":
    comparetoref()

hetsites()





