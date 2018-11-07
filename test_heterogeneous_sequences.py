#!/usr/bin/env python3

from unittest import TestCase
from dark.reads import Read, Reads
from collections import Counter, defaultdict
from heterogeneous_sequences import heterogeneousSites, compareToRef


class TestHeterogeneousSites(TestCase):
    """
    Test the heterogeneousSites function.
    """

    def testOneRead(self):
        read = Read('id', 'ACCG')
        reads = Reads([read])

        self.assertEqual({}, heterogeneousSites(reads, len(read), 1))

    def testHomogeneousReads(self):
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCG')])

        self.assertEqual({}, heterogeneousSites(reads, len(read), 1))

    def testHeterogeneousReadsOneDifference(self):
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCC')])

        self.assertEqual({3: Counter({'G': 1, 'C': 1})}, heterogeneousSites(reads, len(read), 1))

    def testHeterogeneousReadsTwoDifferences(self):
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'TCCC')])

        self.assertEqual({0: Counter({'A': 1, 'T': 1}), 3: Counter({'G': 1, 'C': 1})}, heterogeneousSites(reads, len(read), 1))

    def testHeterogeneousReadsFractionHigh(self):
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCC'), Read('id3', 'ACCC')])

        self.assertEqual({3: Counter({'G': 1, 'C': 2})}, heterogeneousSites(reads, len(read), 0.7))

    def testHeterogeneousReadsFractionLow(self):
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCC'), Read('id3', 'ACCC')])

        self.assertEqual({}, heterogeneousSites(reads, len(read), 0.6))

    def testHeterogeneousReadsFractionLowWithOneDifference(self):
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'TCCG'), Read('id3', 'TCCG'), Read('id4', 'ACCG')])

        self.assertEqual({0: Counter({'A': 2, 'T': 2})}, heterogeneousSites(reads, len(read), 0.6))

    def testGaps(self):
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', '-CC-')])

        self.assertEqual({}, heterogeneousSites(reads, len(read), 1))

    def testEmptyReads(self):
        reads = Reads()

        with self.assertRaisesRegex(NameError, ''):
            heterogeneousSites(reads, len(read), 1)

class TestCompareToReference(TestCase):
    """
    Test the compareToRef function.
    """

    def testEmtpyReference(self):
        referenceread = Read('', '')

        read = Read('id', 'ACCG')
        reads = Reads([read])

        with self.assertRaisesRegex(ValueError, 'Reference read length is not the same as length of other reads.'):
            compareToRef(referenceread, reads, len(read))

    def testNoRead(self):
        referenceread = Read('idref', 'ACCG')

        reads = Reads()

        self.assertEqual({}, compareToRef(referenceread, reads, len(referenceread)))

    def testOneReadNoDifference(self):
        referenceread = Read('idref', 'ACCG')

        read = Read('id', 'ACCG')
        reads = Reads([read])

        self.assertEqual({}, compareToRef(referenceread, reads, len(read)))

    def testOneReadOneDifference(self): 
        referenceread = Read('idref', 'ACCG')

        read = Read('id', 'ACCA')
        reads = Reads([read])

        self.assertEqual({'G to A': 1}, compareToRef(referenceread, reads, len(read)))

    def testAllPossibleDifferences(self):
        referenceread = Read('idref', 'AAACCCGGGTTT')

        read = Read('id', 'CGTAGTACTACG')
        reads = Reads([read])

        self.assertEqual({'A to C': 1, 'A to G': 1, 'A to T': 1, 'C to A': 1, 'C to G': 1, 
            'C to T': 1, 'G to A': 1, 'G to C': 1, 'G to T': 1,'T to A': 1, 'T to C': 1, 
            'T to G': 1}, compareToRef(referenceread, reads, len(read)))

    def testTwoReadsTenDifferences(self): 
        referenceread = Read('idref', 'ACCGGG')

        read = Read('id', 'AGGAAA')
        reads = Reads([read, Read('id2', 'AGGAAA')])

        self.assertEqual({'C to G': 4, 'G to A': 6}, compareToRef(referenceread, reads, len(read)))

    def testGapInReference(self):
        referenceread = Read('idref', 'ACC-')

        read = Read('id', 'ACCG')
        reads = Reads([read])

        self.assertEqual({}, compareToRef(referenceread, reads, len(read)))

    def testGapInReadNoMutation(self):
        referenceread = Read('idref', 'ACCG')

        read = Read('id', 'ACC-')
        reads = Reads([read])

        self.assertEqual({}, compareToRef(referenceread, reads, len(read)))

    def testGapInReadAndMutation(self):
        referenceread = Read('idref', 'ACCGG')

        read = Read('id', 'ACC-G')
        reads = Reads([read, Read('id2', 'ACCAA')])

        self.assertEqual({'G to A': 2}, compareToRef(referenceread, reads, len(read)))

    def testLengthDifference(self):
        referenceread = Read('idref', 'ACC')

        read = Read('id', 'ACCG')
        reads = Reads([read])

        with self.assertRaisesRegex(ValueError, 'Reference read length is not the same as length of other reads.'):
            compareToRef(referenceread, reads, len(read))



