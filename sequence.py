#!/usr/bin/env python

#   Copyright (C) 2023-2024 Simon Crase

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

'''Chapter 2: Which DNA Patterns play the role of molecular clocks'''

from random import randint, random
from sys import float_info
from unittest import main, TestCase, skip
import numpy as np
from rosalind import hamm, k_mers
from reference_tables import bases


def enumerateMotifs(k,d,Dna):
    '''
     BA2A 	Implement MotifEnumeration

    Parameters:
        k        Length of the motifs that we are looking for
        d        Maximum allowable distance
        Dna      The strings to be searched

    Returns:
        All (k, d)-motifs in Dna, i.e. all k-mer that appear in every string from Dna with at most d mismatches.
    '''
    class ApproximateMatches:
        '''
        This class is used to build a dictionary of near matches for kmers
        '''
        def __init__(self):
            self.matches = {}

        def add_matches(self,kmer):
            '''
            Add other kmers to dictionary entry for supplied kmer if within d
            '''
            if not kmer in self.matches:
                self.matches[kmer] = [km for km in kmers if hamm(km,kmer) <=d ]

        def get_matches(self,kmer):
            return self.matches[kmer]

    def close_enough(pattern):
        '''
        Verify that a pattern is within specified maximum distance of every string in Dna
        '''
        def match(string):
            '''
            Verify that a pattern is within specified maximum distance of some substring
            '''
            for i in range(len(string)-k+1):
                if hamm(string[i:i+k],pattern) <= d:
                    return True
            return False

        for string in Dna:
            if not match(string):
                return False
        return True

    approximate_matches = ApproximateMatches()
    kmers = k_mers(k)
    patterns = set()
    for string in Dna:
        for i in range(len(string)-k+1):
            kmer = string[i:i+k]
            approximate_matches.add_matches(kmer)
            for pattern in approximate_matches.get_matches(kmer):
                if close_enough(pattern) and pattern not in patterns:
                    patterns.add(pattern)
                    yield pattern


def medianString(k,dna):
    '''
    BA2B 	Find a Median String

    Input: An integer k and a collection of strings Dna.

    Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers
    '''
    distance = float_info.max
    closest_k_mer = None
    for candidate_k_mer in k_mers(k):
        candidate_distance = distanceBetweenPatternAndStrings(candidate_k_mer,dna)
        if distance > candidate_distance:
            distance = candidate_distance
            closest_k_mer = candidate_k_mer
    return closest_k_mer

def mostProbable(text,n,profile):
    '''
     BA2C Find a Profile-most Probable k-mer in a String

     Input: A string Text, an integer k, and a 4 × k matrix Profile.

     Return: A Profile-most probable k-mer in Text.
    '''
    def log_prob(kmer):
        '''log(probability of kmer given profile)'''
        def log(x):
            '''
            Workaround for #149. This suppresses the warning, but does not affect result of mostProbable,
            since result is computed using argmax
            '''
            return np.log(x) if x > 0 else -np.inf
        return sum([log(profile[bases.find(kmer[j])][j]) for j in range(n)])

    def findMostProbable():
        kmers = [text[i:i+n] for i in range(len(text)-n+1)]
        return kmers[np.argmax( [log_prob(s) for s in kmers])]

    return findMostProbable()


def greedyMotifSearch(k,t,dna,
                      pseudo_counts=False):
    '''
    BA2D 	Implement GreedyMotifSearch
    BA2E 	Implement GreedyMotifSearch with Pseudocounts

    Input: Integers k and t, followed by a collection of strings Dna.
           Optional parameter pseudo_counts specifies whether pseudo counts are to be used

    Return: A collection of strings BestMotifs resulting from running
    GreedyMotifSearch(Dna, k, t). If at any step you find more than one
    Profile-most probable k-mer in a given string, use the one occurring first.
    '''
    def ones_or_zeros(pseudo_counts):
        '''
        Used to initialize counts to 1s or 0s, depending on pseudo_counts
        '''
        return np.ones((len(bases),k),dtype=int) if pseudo_counts else np.zeros((len(bases),k),dtype=int)

    def count_occurrences_of_bases(motifs):
        '''
        Create an array containing the count of occurences of
        each base at each position, summed over all motifs
        '''
        counts = ones_or_zeros(pseudo_counts)
        for kmer in motifs:
            for j in range(k):
                counts[bases.find(kmer[j]),j] += 1
        return counts

    def profile(motifs):
        return count_occurrences_of_bases(motifs)/float(len(motifs))

    def score(motifs):
        counts = count_occurrences_of_bases(motifs)
        total = 0
        for j in range(k):
            total += (len(bases) - counts[:,j].max())
        return total

    bestMotifs = [genome[0:k] for genome in dna]
    for motif in [dna[0][i:i+k] for i in range(len(dna[0])-k+1)]:
        motifs = [motif]
        for i in range(1,t):
            motifs.append(mostProbable(dna[i],k,profile(motifs)))
        if score(motifs) < score(bestMotifs):
            bestMotifs = motifs
    return bestMotifs

def randomized_motif_search(k,t,dna,eps=1):
    '''
    BA2F Implement RandomizedMotifSearch

    Parameters:
        k
        t
        dna
        eps
    '''
    def score(k,motifs):
        total = 0
        for j in range(k):
            counts = np.zeros(len(bases),dtype=np.int32)
            for motif in motifs:
                i = bases.find(motif[j])
                counts[i] += 1
            max_count = -1
            ii = -1
            for i in range(len(bases)):
                if max_count < counts[i]:
                    ii = i
                    max_count = counts[ii]
            for i in range(len(bases)):
                if i != ii:
                    total += counts[i]
        return total

    def counts(motifs):
        matrix = np.ones((len(bases),k),dtype=int)
        for i in range(len(bases)):
            for j in range(k):
                matrix[i,j] *= eps
        for kmer in motifs:
            for j in range(k):
                i = bases.find(kmer[j])
                matrix[i,j] += 1
        return matrix

    def Motifs(profile,dna):
        def get_k(Profile):
            return len(Profile[0])
        def prob(kmer):
            p = 1
            for j in range(k):
                i = bases.find(kmer[j])
                p *= profile[i][j]
            return p
        k = get_k(profile)
        motifs = []
        for s in dna:
            max_probability =- 1
            most_probable_kmer = ''
            for kmer in [s[i:i+k].upper() for i in range(len(s)-k+1)]:
                if max_probability < prob(kmer):
                    max_probability = prob(kmer)
                    most_probable_kmer = kmer
            motifs.append(most_probable_kmer)
        return motifs

    def Profile(motifs):
        matrix = counts(motifs)
        probabilities = np.zeros((len(bases),k),dtype=float)
        for i in range(len(bases)):
            for j in range(k):
                probabilities[i,j] = matrix[i,j]/float(len(motifs))
        return probabilities

    def random_kmer(string):
        i = randint(0,len(string)-k)
        return string[i:i+k]

    motifs = []

    for i in range(t):
        motifs.append(random_kmer(dna[i]))
    bestMotifs = motifs
    while True:
        profile = Profile(motifs)
        motifs = Motifs(profile, dna)
        if score(k,motifs) < score(k,bestMotifs):
            bestMotifs = motifs
        else:
            return (score(k,bestMotifs),bestMotifs)

def randomized_motif_search_driver(k,t,dna,N=2000):
    '''
    BA2F Implement RandomizedMotifSearch

    Repeatedly execute randomized_motif_search and choose best scoring solution
    '''
    best_score = float_info.max
    best_motifs = []
    for i in range(N):
        (score,candidate_motifs) = randomized_motif_search(k,t,dna)
        if score < best_score:
            best_score = score
            best_motifs = candidate_motifs

    return (best_score,best_motifs)


def gibbs(k,t,dna,
          n = 20,
          eps = 1):
    '''
    BA2G 	Implement GibbsSampler

     Input: Integers k, t, and N, followed by a collection of strings Dna.

     Return: The strings BestMotifs resulting from running
     GibbsSampler(Dna, k, t, N) with 20 random starts.
    '''
    def score(k,motifs):
        total = 0
        for j in range(k):
            counts = np.zeros(len(bases),dtype=np.int32)
            for motif in motifs:
                counts[bases.find(motif[j])] += 1
            max_count =- 1
            for i in range(len(bases)):
                if max_count < counts[i]:
                    ii = i
                    max_count = counts[ii]
            # index_max_count = np.argmax(counts)
            # total += sum([counts[i] for i in range(len(bases)) if i != index_max_count])
            for i in range(len(bases)):
                if i != ii:
                    total += counts[i]
        return total

    def random_kmer(string):
        i = randint(0,len(string)-k)
        return string[i:i+k]

    def dropOneMotif(motifs,i):
        return [motifs[j] for j in range(len(motifs)) if j!=i]

    def counts(motifs):
        matrix = np.ones((len(bases),k),dtype=int)
        for i in range(len(bases)):
            for j in range(k):
                matrix[i,j] *= eps
        for kmer in motifs:
            for j in range(k):
                i = bases.find(kmer[j])
                matrix[i,j] += 1
        return matrix

    def Profile(motifs):
        matrix = counts(motifs)
        probabilities = np.zeros((len(bases),k),dtype=float)
        for i in range(len(bases)):
            for j in range(k):
                probabilities[i,j] = matrix[i,j]/float(len(motifs))
        return probabilities

    def probability(kmer,profile):
        p = 1
        for j in range(len(kmer)):
            i = bases.find(kmer[j])
            p *= profile[i][j]
        return p

    def accumulate(probabilities):
        total = 0
        cumulative = []
        for p in probabilities:
            total += p
            cumulative.append(total)
        return cumulative

    def generate(probabilities):
        accumulated=accumulate(probabilities)
        rr = accumulated[len(accumulated)-1]*random()
        i = 0
        while accumulated[i] <= rr:
            i += 1
        return i

    motifs = []

    for i in range(t):
        motifs.append(random_kmer(dna[i]))
    bestMotifs = motifs

    trace = []
    best_score = float_info.max
    for j in range(n):
        i = randint(0,t-1)
        profile = Profile(dropOneMotif(motifs,i))
        motif_index = generate([probability(dna[i][ll:ll+k],profile)
                              for ll in range(len(dna[i])-k)])
        motifs[i] = dna[i][motif_index:motif_index+k]
        sc = score(k,motifs)
        if  sc < best_score:
            best_score = sc
            bestMotifs = motifs
        trace.append(best_score)

    return (score(k,bestMotifs),bestMotifs,trace)


def distanceBetweenPatternAndStrings (pattern,dna):
    '''
    BA2H 	Implement DistanceBetweenPatternAndStrings
    '''

    def hamming(pattern,genome):
        '''Extend Hamming distance to work with string of unequal length'''
        return min([hamm(pattern,genome[i:i+len(pattern)])
                    for i in range(len(genome)-len(pattern)+1)])
    return sum([hamming(pattern,motif) for motif in dna])

if __name__=='__main__':
    class Test_2_Sequence(TestCase):

        def test_ba2a_light(self):
            '''
            BA2A Implement MotifEnumeration

            Test with sample data
            '''
            self.assertEqual(['ATA', 'ATT', 'GTT', 'TTT'],
                             sorted(
                                 list(
                                     enumerateMotifs(3,
                                                 1,
                                                 ['ATTTGGC',
                                                  'TGCCTTA',
                                                  'CGGTATC',
                                                  'GAAAATT']))))

        def test_ba2a_rosalind(self):
            '''
            BA2A Implement MotifEnumeration

            Test with a Rosalind dataset
            '''
            self.assertEqual([
                            'AAAGC', 'AACAG', 'AACTG', 'AAGCG', 'AAGCT', 'AAGGC', 'AAGGG', 'AATAG', 'AATGG', 'ACAAG',
                            'ACACG', 'ACAGG', 'ACAGT', 'ACATG', 'ACATT', 'ACCGG', 'ACCTC', 'ACGAC', 'ACGAT', 'ACGCA',
                            'ACGCC', 'ACGCG', 'ACGCT', 'ACGGC', 'ACGGG', 'ACGGT', 'ACGTG', 'ACGTT', 'ACTCA', 'ACTCC',
                            'ACTCG', 'ACTGA', 'ACTGG', 'ACTGT', 'ACTTC', 'AGAAT', 'AGACC', 'AGAGG', 'AGAGT', 'AGATA',
                            'AGATC', 'AGATT', 'AGCAT', 'AGCCT', 'AGCGA', 'AGCGC', 'AGCGG', 'AGCGT', 'AGCTC', 'AGCTG',
                            'AGGAC', 'AGGAT', 'AGGCC', 'AGGGA', 'AGGGC', 'AGGGT', 'AGGTC', 'AGTCC', 'AGTCT', 'AGTGA',
                            'AGTGG', 'ATAAC', 'ATAAG', 'ATAGA', 'ATAGC', 'ATAGT', 'ATATC', 'ATCAC', 'ATCCA', 'ATCGA',
                            'ATCGC', 'ATCGG', 'ATCGT', 'ATCTC', 'ATCTG', 'ATGAT', 'ATGCC', 'ATGCG', 'ATGCT', 'ATGGC',
                            'ATGGT', 'ATTCA', 'ATTCG', 'ATTCT', 'ATTGA', 'ATTGG', 'CAAAG', 'CAAAT', 'CAACC', 'CAAGA',
                            'CAAGG', 'CAATG', 'CACAA', 'CACAG', 'CACAT', 'CACCG', 'CACCT', 'CACGA', 'CACGG', 'CACTG',
                            'CAGAA', 'CAGAT', 'CAGCA', 'CAGGA', 'CAGGG', 'CAGTA', 'CAGTC', 'CAGTG', 'CAGTT', 'CATAG',
                            'CATCA', 'CATCG', 'CATGA', 'CATGC', 'CATGG', 'CATTC', 'CATTG', 'CCAAT', 'CCACA', 'CCACC',
                            'CCACG', 'CCAGG', 'CCATA', 'CCATC', 'CCATG', 'CCCAA', 'CCCCG', 'CCCGT', 'CCCTA', 'CCCTG',
                            'CCGAA', 'CCGAT', 'CCGCG', 'CCGCT', 'CCGGA', 'CCGGT', 'CCGTC', 'CCGTG', 'CCGTT', 'CCTAA',
                            'CCTAT', 'CCTCA', 'CCTCG', 'CCTGA', 'CCTGT', 'CCTTC', 'CCTTG', 'CGAAC', 'CGAAG', 'CGACA',
                            'CGACC', 'CGACT', 'CGAGA', 'CGAGC', 'CGAGG', 'CGATC', 'CGATT', 'CGCAA', 'CGCAG', 'CGCCA',
                            'CGCCT', 'CGCGA', 'CGCGT', 'CGCTA', 'CGCTG', 'CGCTT', 'CGGAC', 'CGGAT', 'CGGCT', 'CGGGA',
                            'CGGTA', 'CGGTC', 'CGGTG', 'CGGTT', 'CGTCA', 'CGTCG', 'CGTCT', 'CGTGA', 'CGTGC', 'CGTGG',
                            'CGTGT', 'CGTTA', 'CGTTC', 'CGTTG', 'CGTTT', 'CTAAA', 'CTAAT', 'CTACA', 'CTACG', 'CTACT',
                            'CTATC', 'CTCAG', 'CTCAT', 'CTCCA', 'CTCGA', 'CTCGC', 'CTCGG', 'CTCGT', 'CTCTC', 'CTCTG',
                            'CTGAA', 'CTGAC', 'CTGAG', 'CTGAT', 'CTGCC', 'CTGGA', 'CTGGC', 'CTGGT', 'CTGTC', 'CTGTG',
                            'CTTCA', 'CTTCG', 'CTTGC', 'CTTTC', 'GAAGC', 'GAATT', 'GACAT', 'GACCA', 'GACCT', 'GAGAT',
                            'GAGCC', 'GAGCG', 'GAGCT', 'GAGGA', 'GAGGG', 'GAGGT', 'GATAC', 'GATAG', 'GATCA', 'GATCC',
                            'GATCG', 'GATCT', 'GATGC', 'GATGG', 'GATTA', 'GATTC', 'GATTG', 'GCAAG', 'GCACC', 'GCACG',
                            'GCACT', 'GCAGT', 'GCATC', 'GCATG', 'GCATT', 'GCCAT', 'GCCCA', 'GCCCG', 'GCCTA', 'GCCTC',
                            'GCCTG', 'GCCTT', 'GCGAA', 'GCGAT', 'GCGCG', 'GCGCT', 'GCGGT', 'GCGTA', 'GCGTC', 'GCGTT',
                            'GCTAT', 'GCTCA', 'GCTCC', 'GCTCG', 'GCTCT', 'GCTGA', 'GCTGC', 'GCTGG', 'GCTGT', 'GCTTG',
                            'GCTTT', 'GGAAT', 'GGACG', 'GGACT', 'GGAGC', 'GGATC', 'GGATG', 'GGATT', 'GGCAC', 'GGCAT',
                            'GGCCG', 'GGCCT', 'GGCTC', 'GGCTG', 'GGCTT', 'GGGAT', 'GGGCG', 'GGGCT', 'GGTAC', 'GGTAT',
                            'GGTCC', 'GGTCG', 'GGTCT', 'GGTGC', 'GGTGG', 'GGTTG', 'GTAAC', 'GTACC', 'GTACT', 'GTAGT',
                            'GTCAC', 'GTCAT', 'GTCCA', 'GTCCC', 'GTCCT', 'GTCGC', 'GTCTA', 'GTCTC', 'GTGAA', 'GTGAC',
                            'GTGAG', 'GTGAT', 'GTTAC', 'GTTAG', 'GTTCA', 'GTTCG', 'GTTGC', 'GTTGG', 'TAAAG', 'TAACA',
                            'TAACC', 'TAACG', 'TAAGC', 'TAAGG', 'TAATG', 'TACGA', 'TACGG', 'TACGT', 'TACTA', 'TACTC',
                            'TACTT', 'TAGAG', 'TAGAT', 'TAGCA', 'TAGCG', 'TAGCT', 'TAGGA', 'TAGGT', 'TAGTC', 'TAGTG',
                            'TAGTT', 'TATAG', 'TATCC', 'TATCG', 'TATGC', 'TATGG', 'TCAAG', 'TCACA', 'TCACC', 'TCACT',
                            'TCAGA', 'TCAGC', 'TCAGT', 'TCATG', 'TCATT', 'TCCAC', 'TCCAT', 'TCCCG', 'TCCCT', 'TCCGG',
                            'TCCTT', 'TCGAA', 'TCGAG', 'TCGAT', 'TCGCA', 'TCGCG', 'TCGCT', 'TCGGA', 'TCGGG', 'TCGGT',
                            'TCGTA', 'TCGTG', 'TCGTT', 'TCTAG', 'TCTAT', 'TCTCC', 'TCTCG', 'TCTCT', 'TCTGA', 'TCTGG',
                            'TGACA', 'TGACG', 'TGACT', 'TGAGC', 'TGATC', 'TGATG', 'TGCAA', 'TGCAC', 'TGCCC', 'TGCCT',
                            'TGCGT', 'TGCTC', 'TGCTG', 'TGGAT', 'TGGCA', 'TGGCC', 'TGGCG', 'TGGCT', 'TGGGT', 'TGGTC',
                            'TGTCG', 'TGTCT', 'TTACC', 'TTACG', 'TTACT', 'TTAGC', 'TTAGG', 'TTATA', 'TTATG', 'TTCAC',
                            'TTCAG', 'TTGAC', 'TTGCA', 'TTGCC', 'TTGCG', 'TTGCT', 'TTGGC', 'TTGGT', 'TTGTC', 'TTTCG'
                            ],
                             sorted(
                                 list(enumerateMotifs(5,
                                                 2,
                                                 ['CATTAACGTGGTCCTCCGATGTAAT',
                                                  'CCTCATGATGTCTTTATTGAACGTT',
                                                  'CCGCTCAGGTCTAGCTTTGGGATAT',
                                                  'TAGAGTCTTAATCGTTCACGGCGGT',
                                                  'GCGCGCTAGCGCATACTGCCCCGGT',
                                                  'CCGATAGGCACTCGCTTAGGCAATA',
                                                  'CACAGGGTTCTCACGCCCCAGCGAT',
                                                  'CAAGTCACGATCGATGACTATCAAT',
                                                  'ATACAGCGGTCGGTGGCCGGATGGC',
                                                  'TTCGACTCCACTAGCGCGATGCGCC']))))


        def test_ba2b_light(self):
            '''BA2B Find a Median String'''
            self.assertEqual('GAC',
                             medianString(3,
                                          [
                                            'AAATTGACGCAT',
                                            'GACGACCACGTT',
                                            'CGTCAGCGCCTG',
                                            'GCTGAGCACCGG',
                                            'AGTACGGGACAG'
                             ]))


        def test_ba2c_light(self):
            '''BA2C Find a Profile-most Probable k-mer in a String'''
            self.assertEqual(
                'CCGAG',
                mostProbable(
                    'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT',\
                    5,
                    [[0.2, 0.2, 0.3, 0.2, 0.3],
                     [0.4, 0.3, 0.1, 0.5, 0.1],
                     [0.3, 0.3, 0.5, 0.2, 0.4],
                     [0.1, 0.2, 0.1, 0.1, 0.2]]))

        def test_ba2d_light(self):
            ''' BA2D Implement GreedyMotifSearch'''
            motifs = greedyMotifSearch(3,
                                       5,
                                       [
                                        'GGCGTTCAGGCA',
                                        'AAGAATCAGTCA',
                                        'CAAGGAGTTCGC',
                                        'CACGTCAATCAC',
                                        'CAATAATATTCG'
                                       ])
            self.assertEqual(['CAA','CAA','CAA','CAG','CAG'],sorted(motifs))


        def test_ba2e_light(self):
            '''BA2E Implement GreedyMotifSearch with Pseudocounts'''
            motifs = greedyMotifSearch(3,
                                       5,
                                       [
                                        'GGCGTTCAGGCA',
                                        'AAGAATCAGTCA',
                                        'CAAGGAGTTCGC',
                                        'CACGTCAATCAC',
                                        'CAATAATATTCG'
                                        ],
                                       pseudo_counts=True)
            self.assertEqual(['ATC','ATC','TTC','TTC','TTC'],sorted(motifs))

        def test_ba2f(self):
            '''
            BA2F Implement RandomizedMotifSearch

            NB: this test fails sometimes (randomness!)
            '''
            (c,x) = randomized_motif_search_driver(8, 5,[
                'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'],1000)

            self.assertIn('TCTCGGGG',x)
            self.assertIn('CCAAGGTG',x)
            self.assertIn('TACAGGCG',x)
            self.assertIn('TTCAGGTG',x)
            self.assertIn('TCCACGTG',x)

        def test_ba2g(self):
            '''
            BA2G Implement GibbsSampler

            We will run Gibbs multiple times, as it does not always find the optimum solution.
            '''
            s0 = float_info.max
            for i in range(1000):
                score,motifs,trace = gibbs(8, 5,
                                               [
                                                'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                                                'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                                                'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                                                'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                                                'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'],
                                               n = 20
                                               )
                if score < s0:
                    s0=score
                    bestMotifs = motifs
                    print (score,bestMotifs)
            self.assertEqual(5,len(bestMotifs))
            self.assertIn('TCTCGGGG', bestMotifs)
            self.assertIn('CCAAGGTG', bestMotifs)
            self.assertIn('TACAGGCG', bestMotifs)
            self.assertIn('TTCAGGTG', bestMotifs)
            self.assertIn('TCCACGTG', bestMotifs)



        def test_ba2h_light(self):
            '''BA2H Implement DistanceBetweenPatternAndStrings'''
            self.assertEqual(
                5,
                distanceBetweenPatternAndStrings('AAA',
                                                 ['TTACCTTAAC',
                                                  'GATATCTGTC',
                                                  'ACGGCGTTCG',
                                                  'CCCTAAAGAG',
                                                  'CGTCAGAGGT']))


    main()
