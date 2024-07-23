#!/usr/bin/env python

#   Copyright (C) 2017-2024 Simon Crase

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

''' Chapter 1: Where in the Genome does DNA replication begin?'''

from unittest import main, skip, TestCase
import numpy as np
from numpy.testing import assert_array_equal
from rosalind import subs, hamm, create_frequency_table, k_mers
from reference_tables import BASES


def countOccurrences(pattern,string):
    '''
    BA1A	Compute the Number of Times a Pattern Appears in a Text

     We define Count(Text, Pattern) as the number of times that a k-mer Pattern
     appears as a substring of Text. For example,

     Count(ACAACTATGCATACTATCGGGAACTATCCT,ACTAT)=3.

     note that Count(CGATATATCCATAG, ATA) is equal to 3 (not 2) since we should
     account for overlapping occurrences of Pattern in Text.

     Input: {DNA strings}} Text and Pattern.

     Return: Count(Text, Pattern).
    '''
    return len(findOccurences(pattern,string))


def find_most_frequent_words(string,k):
    '''
    BA1B	Find the Most Frequent Words in a String

    We say that Pattern is a most frequent k-mer in Text if it maximizes
    Count(Text, Pattern) among all k-mers. For example, "ACTAT" is a most
    frequent 5-mer in "ACAACTATGCATCACTATCGGGAACTATCCT", and "ATA" is a most
    frequent 3-mer of "CGATATATCCATAG".

    Input: A DNA string Text and an integer k.

    Output: All most frequent k-mers in Text (in any order).
    '''
    most_frequent_words = []
    max_k = -1
    for kmer,count in create_frequency_table(string,k).items():
        if count>max_k:
            max_k = count
            most_frequent_words = []
        if count == max_k:
            most_frequent_words.append(kmer)
    return most_frequent_words


def findOccurences(pattern,string):
    '''
    BA1D	Find All Occurrences of a Pattern in a String

    Input: Strings Pattern and Genome.

    Return: All starting positions in Genome where Pattern appears as a substring.
            Use 0-based indexing.
    '''
    return [pos-1 for pos in subs(string,pattern)]




def findClumps(genome,k,L,t):
    '''
     BA1E	Find Patterns Forming Clumps in a String

    Given integers L and t, a string Pattern forms an (L, t)-clump inside a
    (larger) string Genome if there is an interval of Genome of length L in which
    Pattern appears at least t times. For example, TGCA forms a (25,3)-clump in
    the following Genome: gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

    Input: A string Genome, and integers k, L, and t.

    Return: All distinct k-mers forming (L, t)-clumps in Genome.
    '''
    def update_patterns(frequencies):
        for (kmer,count) in frequencies.items():
            if count>=t:
                patterns.append(kmer)

    patterns = []
    frequencies = create_frequency_table(genome[0:L],k)
    update_patterns(frequencies)
    for i in range(1,len(genome)-L+1):
        head = genome[i-1:i+k-1]
        frequencies[head] -= 1
        tail = genome[i+L-k:i+L]
        if tail in frequencies:
            frequencies[tail] += 1
        else:
            frequencies[tail] = 1
        update_patterns(frequencies)

    return list(set(patterns))



def find_minimum_skew(genome):
    '''
    # BA1F	Find a Position in a Genome Minimizing the Skew

    Define the skew of a DNA string Genome, denoted Skew(Genome), as the
    difference between the total number of occurrences of 'G' and 'C' in Genome.

    Parameters:
        genome A DNA string

    Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of
            i (from 0 to |Genome|).
    '''
    SKEW_STEP={
        'A':0,
        'C':-1,
        'G': +1,
        'T': 0
    }
    positions = []
    min_skew = 2
    skew = 0
    pos = 0
    for nucleotide in genome:
        pos += 1
        skew += SKEW_STEP[nucleotide]
        if min_skew > skew:
            min_skew = skew
            positions = [pos]
        elif min_skew == skew:
            positions.append(pos)

    return positions, min_skew


def findApproximateOccurrences(pattern,text,d):
    '''
     BA1H	Find All Approximate Occurrences of a Pattern in a String

    Input: Strings Pattern and Text along with an integer d.

    Return: All starting positions where Pattern appears as a substring of Text
            with at most d mismatches.

    '''
    return [i
            for i in range(len(text)-len(pattern)+1)
            if hamm(pattern,text[i:i+len(pattern)])<=d]



def find_mismatches(pattern,text,d):
    ''' helper for BA1I'''
    return findApproximateOccurrences(pattern,text,d)

def find_mismatches_and_rc(pattern,text,d):
    '''helper for BA1J'''
    return findApproximateOccurrences(pattern,text,d) + findApproximateOccurrences(revc(pattern),text,d)



def findMostFrequentWordsWithMismatches(text,k,d,
                                        find=find_mismatches):
    '''
     BA1I	Find the Most Frequent Words with Mismatches in a String
     BA1J	Find Frequent Words with Mismatches and Reverse Complements

     A most frequent k-mer with up to d mismatches in Text is simply a string
     Pattern maximizing Countd(Text, Pattern) among all k-mers. Note that Pattern
     does not need to actually appear as a substring of Text; for example, AAAAA
     is the most frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG,
     even though AAAAA does not appear exactly in this string.

     Input: A string Text as well as integers k and d.

     Return: All most frequent k-mers with up to d mismatches in Text.
    '''
    matches = []
    max_count = -1
    for pattern in k_mers(k):
        count = len(find(pattern,text,d))
        if count > max_count:
            max_count = count
            matches = []
        if count == max_count:
            matches.append(pattern)

    return (max_count,matches)

def generateFrequencyArray(text,k):
    '''
    BA1K	Generate the Frequency Array of a Strings
    Given an integer k, we define the frequency array of a string Text as an
    array of length 4**k, where the i-th element of the array holds the number of
    times that the i-th k-mer (in the lexicographic order) appears in Text.

    Input: A DNA string Text and an integer k.

    Return: The frequency array of k-mers in Text.
    '''
    frequencies = np.zeros((4**k))
    for i in range(len(text)-k+1):
        frequencies[patternToNumber(text[i:i+k])] += 1
    return frequencies


def patternToNumber(kmer):
    '''BA1L	Implement PatternToNumber'''
    n = 0
    for letter in kmer:
        n *= 4
        n += BASES.find(letter)
    return n

def numberToPattern(n,k):
    '''BA1M	Implement NumberToPattern'''
    pattern = ''
    nn = n
    for i in range(k):
        pattern = pattern+bases[nn%4]
        nn //= 4
    return pattern[::-1]

def generate_dNeighborhood(pattern,d):
    '''
    BA1N Generate the d-Neighborhood of a String

    The d-neighborhood Neighbors(Pattern, d) is the set of all k-mers whose
    Hamming distance from Pattern does not exceed d.

    Input: A DNA string Pattern and an integer d.

    Return: The collection of strings Neighbors(Pattern, d).
    '''
    def neighbours(p):
        neighbours = []
        for i in range(len(p)):
            for base in BASES:
                neighbours.append(p[0:i] + base + p[i+1:])
        return neighbours
    if d==0:
        return [pattern]
    else:
        neighbourhood = []
        start = generate_dNeighborhood(pattern,d-1)
        for p in start:
            for n in neighbours(p):
                neighbourhood.append(n)
        return list(set(start + neighbourhood))

if __name__=='__main__':
    class TestCase_1_Replication(TestCase):
        ''' Tests from Chapter 1: Where in the Genome does DNA replication begin?'''
        def test_ba1a(self):
            '''BA1A	Compute the Number of Times a Pattern Appears in a Text'''
            self.assertEqual(2,countOccurrences('GCG','GCGCG'))

        def test_ba1b1(self):
            '''BA1B	Find the Most Frequent Words in a String'''
            most_frequent_words=find_most_frequent_words('ACGTTGCATGTCGCATGATGCATGAGAGCT',4)
            self.assertIn('CATG',most_frequent_words)
            self.assertIn('GCAT',most_frequent_words)

        def test_ba1b2(self):
            '''BA1B	Find the Most Frequent Words in a String'''
            most_frequent_words= find_most_frequent_words('ACATCTGGCGCCCATCGCCCATCTGACGGTTCTGACGGTTCTGACGGTTAGCAGCTCTGACGGTTCTGACGGTTCTGACGGTTCTGACGGTTCTGACGGTTTCGCGACCGCCCATCTGACGGTTCGCCCATACATCTGGTCGCGACCTGACGGTTACATCTGGCGCCCATCTGACGGTTCGCCCATAGCAGCTCTGACGGTTTCGCGACCGCCCATACATCTGGCTGACGGTTTCGCGACTCGCGACACATCTGGAGCAGCTACATCTGGTCGCGACCTGACGGTTAGCAGCTCGCCCATAGCAGCTCTGACGGTTAGCAGCTCTGACGGTTACATCTGGACATCTGGCGCCCATCGCCCATTCGCGACACATCTGGAGCAGCTCGCCCATTCGCGACTCGCGACAGCAGCTACATCTGGTCGCGACACATCTGGCGCCCATCTGACGGTTACATCTGGTCGCGACACATCTGGCGCCCATTCGCGACACATCTGGACATCTGGAGCAGCTACATCTGGTCGCGACAGCAGCTACATCTGGTCGCGACCTGACGGTTACATCTGGCTGACGGTTCGCCCATCGCCCATCGCCCATTCGCGACAGCAGCTTCGCGACCTGACGGTTACATCTGGCGCCCATACATCTGGTCGCGACAGCAGCTTCGCGACTCGCGACCTGACGGTTTCGCGACACATCTGGCTGACGGTTCTGACGGTTCGCCCATAGCAGCTCTGACGGTTCTGACGGTTAGCAGCTACATCTGGAGCAGCTCGCCCATCGCCCATAGCAGCTAGCAGCTTCGCGACACATCTGGCGCCCATCTGACGGTTTCGCGACAGCAGCT',14)
            self.assertIn('CGGTTCTGACGGTT',most_frequent_words)
            self.assertIn('GACGGTTCTGACGG',most_frequent_words)
            self.assertIn('TGACGGTTCTGACG',most_frequent_words)
            self.assertIn('ACGGTTCTGACGGT',most_frequent_words)
            self.assertIn('CTGACGGTTCTGAC',most_frequent_words)

        def test_ba1b3(self):
            '''BA1B	Find the Most Frequent Words in a String'''
            most_frequent_words=find_most_frequent_words('CGGAAGCGAGATTCGCGTGGCGTGATTCCGGCGGGCGTGGAGAAGCGAGATTCATTCAAGCCGGGAGGCGTGGCGTGGCGTGGCGTGCGGATTCAAGCCGGCGGGCGTGATTCGAGCGGCGGATTCGAGATTCCGGGCGTGCGGGCGTGAAGCGCGTGGAGGAGGCGTGGCGTGCGGGAGGAGAAGCGAGAAGCCGGATTCAAGCAAGCATTCCGGCGGGAGATTCGCGTGGAGGCGTGGAGGCGTGGAGGCGTGCGGCGGGAGATTCAAGCCGGATTCGCGTGGAGAAGCGAGAAGCGCGTGCGGAAGCGAGGAGGAGAAGCATTCGCGTGATTCCGGGAGATTCAAGCATTCGCGTGCGGCGGGAGATTCAAGCGAGGAGGCGTGAAGCAAGCAAGCAAGCGCGTGGCGTGCGGCGGGAGAAGCAAGCGCGTGATTCGAGCGGGCGTGCGGAAGCGAGCGG',12)
            self.assertIn('CGGCGGGAGATT',most_frequent_words)
            self.assertIn('CGGGAGATTCAA',most_frequent_words)
            self.assertIn('CGTGCGGCGGGA',most_frequent_words)
            self.assertIn('CGTGCGGCGGGA',most_frequent_words)
            self.assertIn('CGTGGAGGCGTG',most_frequent_words)
            self.assertIn('CGTGGCGTGCGG',most_frequent_words)
            self.assertIn('GCGTGCGGCGGG',most_frequent_words)
            self.assertIn('GCGTGGAGGCGT',most_frequent_words)
            self.assertIn('GCGTGGCGTGCG',most_frequent_words)
            self.assertIn('GGAGAAGCGAGA',most_frequent_words)
            self.assertIn('GGAGATTCAAGC',most_frequent_words)
            self.assertIn('GGCGGGAGATTC',most_frequent_words)
            self.assertIn('GGGAGATTCAAG',most_frequent_words)
            self.assertIn('GTGCGGCGGGAG',most_frequent_words)
            self.assertIn('TGCGGCGGGAGA',most_frequent_words)

        def test_ba1e(self):
            clumps=findClumps('CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAAC'
                              'ACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC',5,75,4)
            self.assertIn('CGACA',clumps)
            self.assertIn('GAAGA',clumps)
            self.assertIn('AATGT',clumps)

        def test_ba1f(self):
            positions,_=find_minimum_skew(
                'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCT'\
                'TGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG')
            self.assertIn(53,positions)
            self.assertIn(97,positions)

        def test_ba1g(self):
            self.assertEqual(7,hamm('GAGCCTACTAACGGGAT','CATCGTAATGACGGCCT'))

        def test_ba1h(self):
            self.assertEqual([6, 7, 26, 27, 78],
                             findApproximateOccurrences('ATTCTGGA',
                                                        'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC',\
                                                           3))

        def test_ba1i(self):
            _,matches=findMostFrequentWordsWithMismatches('ACGTTGCATGTCGCATGATGCATGAGAGCT',4,1)
            self.assertIn('GATG',matches)
            self.assertIn('ATGC',matches)
            self.assertIn('ATGT',matches)

        def test_ba1k(self):
            assert_array_equal(np.array([2,1,0,0,0,0,2,2,1,2,1,0,0,1,1,0]),
                             generateFrequencyArray('ACGCGGCTCTGAAA',2))

        def test_ba1n(self):
            neighbours=generate_dNeighborhood('ACG',1)
            self.assertIn('CCG', neighbours)
            self.assertIn('TCG', neighbours)
            self.assertIn('GCG', neighbours)
            self.assertIn('AAG', neighbours)
            self.assertIn('ATG', neighbours)
            self.assertIn('AGG', neighbours)
            self.assertIn('ACA', neighbours)
            self.assertIn('ACC', neighbours)
            self.assertIn('ACT', neighbours)
            self.assertIn('ACG', neighbours)

    main()
