#!/usr/bin/env python3

from unittest import TestCase
from dark.reads import Read, Reads
from heterogeneous_sequences import heterogeneousSites, compareToRef


class TestHeterogeneousSites(TestCase):
    """
    Test the heterogeneousSites function.
    """

    def testOneRead(self):
        """
        heterogeneousSites must return an empty dictionary if only one read is given.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read])

        self.assertEqual({}, heterogeneousSites(reads, len(read), 1))

    def testHomogeneousReads(self):
        """
        heterogeneousSites must return an empty dictionary if homogenous reads are given.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCG')])

        self.assertEqual({}, heterogeneousSites(reads, len(read), 1))

    def testHeterogeneousReadsOneDifference(self):
        """
        heterogeneousSites must return a dictionary with one entry as expected if 
        reads given differ at one site.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCC')])

        self.assertEqual({3: {'G': 1, 'C': 1}}, heterogeneousSites(reads, len(read), 1))

    def testHeterogeneousReadsTwoDifferences(self):
        """
        heterogeneousSites must return a dictionary with two entries as expected if 
        reads given differ at two sites.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'TCCC')])

        self.assertEqual({0: {'A': 1, 'T': 1}, 3: {'G': 1, 'C': 1}}, heterogeneousSites(reads, len(read), 1))

    def testHeterogeneousReadsFractionHigh(self):
        """
        heterogeneousSites must return a dictionary with one entry as expected if 
        reads given differ and are less homogeneous than specified by the homogeneity 
        cutoff fraction.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCC'), Read('id3', 'ACCC')])

        self.assertEqual({3: {'G': 1, 'C': 2}}, heterogeneousSites(reads, len(read), 0.7))

    def testHeterogeneousReadsFractionLow(self):
        """
        heterogeneousSites must return an empty dictionary if reads given differ and 
        are more homogeneous than specified by the homogeneity cutoff fraction.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCC'), Read('id3', 'ACCC')])

        self.assertEqual({}, heterogeneousSites(reads, len(read), 0.6))

    def testHeterogeneousReadsFractionLowWithOneDifference(self):
        """
        heterogeneousSites must return a dictionary with one entry if reads given differ 
        at two sites and at one site are more homogeneous than specified by the homogeneity 
        cutoff fraction; at the other site less homogeneous than specified by the 
        homogeneity cutoff fraction.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'TCCG'), Read('id3', 'TCCG'), Read('id4', 'ACCG')])

        self.assertEqual({0: {'A': 2, 'T': 2}}, heterogeneousSites(reads, len(read), 0.6))

    def testGaps(self):
        """
        heterogeneousSites must return an empty dictionary if reads given differ only by gaps;
        gaps do not count towards heterogeneity.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', '-CC-')])

        self.assertEqual({}, heterogeneousSites(reads, len(read), 1))

    def testEmptyReads(self):
        """
        heterogeneousSites must return an empty dictionary if no reads are given.
        """
        reads = Reads()

        self.assertEqual({}, heterogeneousSites(reads, 0, 1))

class TestCompareToReference(TestCase):
    """
    Test the compareToRef function.
    """

    def testEmtpyReference(self):
        """
        compareToRef must raise a ValueError if an empty reference is given.
        """
        referenceread = Read('', '')

        read = Read('id', 'ACCG')
        reads = Reads([read])

        with self.assertRaisesRegex(ValueError, 'Reference read length is not the same as length of other reads.'):
            compareToRef(referenceread, reads, len(read))

    def testNoRead(self):
        """
        compareToRef must return an empty dictionary if no reads are given.
        """
        referenceread = Read('idref', 'ACCG')

        reads = Reads()

        self.assertEqual({'A to C': 0, 
                        'A to G': 0, 
                        'A to T': 0, 
                        'C to A': 0, 
                        'C to G': 0, 
                        'C to T': 0, 
                        'G to A': 0, 
                        'G to C': 0, 
                        'G to T': 0,
                        'T to A': 0, 
                        'T to C': 0, 
                        'T to G': 0}, compareToRef(referenceread, reads, len(referenceread)))

    def testOneReadNoDifference(self):
        """
        compareToRef must return an empty dictionary if reference and read match.
        """
        referenceread = Read('idref', 'ACCG')

        read = Read('id', 'ACCG')
        reads = Reads([read])

        self.assertEqual({'A to C': 0, 
                        'A to G': 0, 
                        'A to T': 0, 
                        'C to A': 0, 
                        'C to G': 0, 
                        'C to T': 0, 
                        'G to A': 0, 
                        'G to C': 0, 
                        'G to T': 0,
                        'T to A': 0, 
                        'T to C': 0, 
                        'T to G': 0}, compareToRef(referenceread, reads, len(read)))

    def testOneReadOneDifference(self): 
        """
        compareToRef must return a dictionary with one entry as expected if read given 
        differs from reference in one instance.
        """
        referenceread = Read('idref', 'ACCG')

        read = Read('id', 'ACCA')
        reads = Reads([read])

        self.assertEqual({'A to C': 0, 
                        'A to G': 0, 
                        'A to T': 0, 
                        'C to A': 0, 
                        'C to G': 0, 
                        'C to T': 0, 
                        'G to A': 1, 
                        'G to C': 0, 
                        'G to T': 0,
                        'T to A': 0, 
                        'T to C': 0, 
                        'T to G': 0}, compareToRef(referenceread, reads, len(read)))

    def testAllPossibleDifferences(self):
        """
        compareToRef must return a dictionary with an entry for each mutationType (12 total)
        if read given differs from reference in all possible ways.
        """
        referenceread = Read('idref', 'AAACCCGGGTTT')

        read = Read('id', 'CGTAGTACTACG')
        reads = Reads([read])

        self.assertEqual({'A to C': 1, 
                        'A to G': 1, 
                        'A to T': 1, 
                        'C to A': 1, 
                        'C to G': 1, 
                        'C to T': 1, 
                        'G to A': 1, 
                        'G to C': 1, 
                        'G to T': 1,
                        'T to A': 1, 
                        'T to C': 1, 
                        'T to G': 1}, 
                        compareToRef(referenceread, reads, len(read)))

    def testTwoReadsTenDifferences(self): 
        """
        compareToRef must return a dictionary with 2 entries as expected if reads given
        differ from the reference in two ways several times.
        """
        referenceread = Read('idref', 'ACCGGG')

        read = Read('id', 'AGGAAA')
        reads = Reads([read, Read('id2', 'AGGAAA')])

        self.assertEqual({'A to C': 0, 
                        'A to G': 0, 
                        'A to T': 0, 
                        'C to A': 0, 
                        'C to G': 4, 
                        'C to T': 0, 
                        'G to A': 6, 
                        'G to C': 0, 
                        'G to T': 0,
                        'T to A': 0, 
                        'T to C': 0, 
                        'T to G': 0}, 
                        compareToRef(referenceread, reads, len(read)))

    def testGapInReference(self):
        """
        compareToRef must return an empty dictionary if reference and read only differ by gaps.
        """
        referenceread = Read('idref', 'ACC-')

        read = Read('id', 'ACCG')
        reads = Reads([read])

        self.assertEqual({'A to C': 0, 
                        'A to G': 0, 
                        'A to T': 0, 
                        'C to A': 0, 
                        'C to G': 0, 
                        'C to T': 0, 
                        'G to A': 0, 
                        'G to C': 0, 
                        'G to T': 0,
                        'T to A': 0, 
                        'T to C': 0, 
                        'T to G': 0}, compareToRef(referenceread, reads, len(read)))

    def testGapInReadNoMutation(self):
        """
        compareToRef must return an empty dictionary if reference and read only differ by gaps.
        """
        referenceread = Read('idref', 'ACCG')

        read = Read('id', 'ACC-')
        reads = Reads([read])

        self.assertEqual({'A to C': 0, 
                        'A to G': 0, 
                        'A to T': 0, 
                        'C to A': 0, 
                        'C to G': 0, 
                        'C to T': 0, 
                        'G to A': 0, 
                        'G to C': 0, 
                        'G to T': 0,
                        'T to A': 0, 
                        'T to C': 0, 
                        'T to G': 0}, compareToRef(referenceread, reads, len(read)))

    def testGapInReadAndMutation(self):
        """
        compareToRef must return a dictionary with one entry if reference and read differ only by gaps 
        at one site, but differ by base at another site.
        """
        referenceread = Read('idref', 'ACCGG')

        read = Read('id', 'ACC-G')
        reads = Reads([read, Read('id2', 'ACCAA')])

        self.assertEqual({'A to C': 0, 
                        'A to G': 0, 
                        'A to T': 0, 
                        'C to A': 0, 
                        'C to G': 0, 
                        'C to T': 0, 
                        'G to A': 2, 
                        'G to C': 0, 
                        'G to T': 0,
                        'T to A': 0, 
                        'T to C': 0, 
                        'T to G': 0}, compareToRef(referenceread, reads, len(read)))

    def testLengthDifference(self):
        """
        compareToRef must reaise a ValueError if the length of the reference and reads doesn't match.
        """
        referenceread = Read('idref', 'ACC')

        read = Read('id', 'ACCG')
        reads = Reads([read])

        with self.assertRaisesRegex(ValueError, 'Reference read length is not the same as length of other reads.'):
            compareToRef(referenceread, reads, len(read))



