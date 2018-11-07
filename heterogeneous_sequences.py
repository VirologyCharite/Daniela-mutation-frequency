#!/usr/bin/env python3

from __future__ import print_function, division
from collections import Counter, defaultdict
from dark.fasta import FastaReads

def heterogeneousSites(reads, length, homogeneous_frequency):
    """
    Returns the sites of heterogeneity in a list of FastaReads of the same length.

    @param reads: An iterable of C{dark.reads.Read} instances.
    @param length: The length of reads.
    @param homogeneous_frequency: Cutoff frequency at which a site is considered to be homogeneous.
    @raises NameError: If C{reads} is empty.
    @return: A C{dict} with C{int} index keys and C{Counter} instances
        as values.
    """
    result = {}

    for index in range(length):
        counts = Counter()

        for read in reads:
            if read.sequence[index] != '-':
                # Ignore gaps.
                counts[read.sequence[index]] += 1

        if counts:
            # The site is not homogenous = more than one nucleotide present.
            if (counts.most_common(1)[0][1] / sum(counts.values()) < homogeneous_frequency):
                result[index] = counts

    return result

def compareToRef(referenceread, reads, length):
    """
    Compares sequences to a reference sequence. Returns a count of mutation types.

    @param referencefasta: A C{FastaRead} instance.
    @param reads: A list of C{FastaReads}
    @param length: The length of reads
    @raises ValueError: if len(referenceread) is not equal to length.
    @return: A C{defaultdict} with C{str} index keys and C{int} instances
        as values.
    """

    reference = referenceread
    totalsitesmutated = 0

    if len(reference.sequence) != length:
        raise ValueError('Reference read length is not the same as length of other reads.')

    mutations = defaultdict(int)

    for index in range(length):
        refbase = reference.sequence[index]
        if refbase == 'U':
            refbase = 'T'       
        elif refbase == '-':
            # Exclude gaps.
            continue

        for read in reads:
            base = read.sequence[index]
            if base == 'U':
                base = 'T'
            if base == '-' or base == refbase:
                # Exclude gaps and non-mutations.
                continue

            mutations['%s to %s' % (refbase, base)] += 1

    return mutations

