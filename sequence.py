#!/usr/bin/env python

#   Copyright (C) 2023 Simon Crase

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

from random   import randint, random
from sys      import float_info
from unittest import main, TestCase, skip

import numpy  as np

from rosalind         import hamm, k_mers
from reference_tables import bases


def enumerateMotifs(k,d,dna):
    '''
     BA2A 	Implement MotifEnumeration

    Input: Integers k and d, followed by a collection of strings Dna.

    Return: All (k, d)-motifs in Dna.
    '''
    approximate_matches = {}
    kmers = k_mers(k)

    def near_matches(kmer):
        if not kmer in approximate_matches:
            approximate_matches[kmer]=[kk for kk in kmers if hamm(kk,kmer)<=d]
        return approximate_matches[kmer]

    def good_enough(pattern):
        def match(string):
            for i in range(len(string)-k+1):
                kmer = string[i:i+k]
                if hamm(kmer,pattern) <= d:
                    return True
            return False

        for string in dna:
            if not match(string):
                return False
        return True

    patterns = []
    for string in dna:
        for i in range(len(string)-k+1):
            kmer = string[i:i+k]
            for pattern in near_matches(kmer):
                if good_enough(pattern):
                    patterns.append(pattern)
    return ' '.join(sorted(list(set(patterns))))


def medianString(k,dna):
    '''
    BA2B 	Find a Median String

    Input: An integer k and a collection of strings Dna.

    Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers
    Pattern. (If multiple answers exist, you may return any one.)

    '''
    def findClosest(d):
        distance = float_info.max
        closest = None
        for k_mer in k_mers(k):
            if distance > d(k_mer,dna):
                distance = d(k_mer,dna)
                closest = k_mer
        return closest
    return findClosest(distanceBetweenPatternAndStrings)




def mostProbable(text,n,profile):
    '''
     BA2C 	Find a Profile-most Probable k-mer in a String

     Input: A string Text, an integer k, and a 4 × k matrix Profile.

     Return: A Profile-most probable k-mer in Text.
    '''
    def log_prob(kmer):
        '''log(probability of kmer given profile)'''
        return sum([np.log(profile[bases.find(kmer[j])][j]) for j in range(n)])

    def findMostProbable():
        kmers = [text[i:i+n] for i in range(len(text)-n+1)]
        return kmers[np.argmax( [log_prob(s) for s in kmers])]

    return findMostProbable()


def greedyMotifSearch(k,t,dna,pseudo_counts=False):
    '''
    BA2D 	Implement GreedyMotifSearch
    BA2E 	Implement GreedyMotifSearch with Pseudocounts

    Input: Integers k and t, followed by a collection of strings Dna.
           Optional parameter pseudo_counts specifies whether pseudo counts are to be used

    Return: A collection of strings BestMotifs resulting from running
    GreedyMotifSearch(Dna, k, t). If at any step you find more than one
    Profile-most probable k-mer in a given string, use the one occurring first.
    '''


    def count_occurrences_of_bases(motifs):
        '''
        Create an array containing the count of occurences of
        each base at each position, summed over all motifs
        '''
        matrix = np.ones((len(bases),k),dtype=int) if pseudo_counts else np.zeros((len(bases),k),dtype=int)
        for kmer in motifs:
            for j in range(k):
                i = bases.find(kmer[j])
                matrix[i,j] += 1
        return matrix

    def profile(motifs):
        return count_occurrences_of_bases(motifs)/float(len(motifs))

    def score(motifs):
        matrix=count_occurrences_of_bases(motifs)
        total=0
        for j in range(k):
            m=0
            for i in range(len(bases)):
                if m<matrix[i,j]:
                    m=matrix[i,j]
            total+=(len(bases)-m)
        return total

    bestMotifs=[genome[0:k] for genome in dna]
    for motif in [dna[0][i:i+k] for i in range(len(dna[0])-k+1)]:
        motifs=[motif]
        for i in range(1,t):
            motifs.append(mostProbable(dna[i],k,profile(motifs)))
        if score(motifs)<score(bestMotifs):
            bestMotifs=motifs
    return bestMotifs





def randomized_motif_search(k,t,dna,eps=1):
    '''
    BA2F 	Implement RandomizedMotifSearch

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

def randomized_motif_search_driver(k,t,dna,N=1000):
    '''
    BA2F Implement RandomizedMotifSearch
    '''
    best = float_info.max
    mm = []
    for i in range(N):
        (sc,motifs) =randomized_motif_search(k,t,dna)
        if sc<best:
            best = sc
            mm = motifs

    return (best,mm)


def gibbs(k,t,n,dna,eps=1):
    '''
    BA2G 	Implement GibbsSampler

     Input: Integers k, t, and N, followed by a collection of strings Dna.

     Return: The strings BestMotifs resulting from running
     GibbsSampler(Dna, k, t, N) with 20 random starts.
     Remember to use pseudocounts
    '''
    def score(k,motifs):
        total=0
        for j in range(k):
            counts=np.zeros(len(bases),dtype=np.int32)
            for motif in motifs:
                i=bases.find(motif[j])
                counts[i]+=1
            max=-1
            ii=-1
            for i in range(len(bases)):
                if max<counts[i]:
                    ii=i
                    max=counts[ii]
            for i in range(len(bases)):
                if i!=ii:
                    total+=counts[i]
        return total

    def random_kmer(string):
        i=randint(0,len(string)-k)
        return string[i:i+k]

    def dropOneMotif(motifs,i):
        return [motifs[j] for j in range(len(motifs)) if j!=i]

    def counts(motifs):
        matrix=np.ones((len(bases),k),dtype=int)
        for i in range(len(bases)):
            for j in range(k):
                matrix[i,j]*=eps
        for kmer in motifs:
            for j in range(k):
                i=bases.find(kmer[j])
                matrix[i,j]+=1
        return matrix

    def Profile(motifs):
        matrix=counts(motifs)
        probabilities=np.zeros((len(bases),k),dtype=float)
        for i in range(len(bases)):
            for j in range(k):
                probabilities[i,j]=matrix[i,j]/float(len(motifs))
        return probabilities

    def probability(kmer,profile):
        p=1
        for j in range(len(kmer)):
            i=bases.find(kmer[j])
            p*=profile[i][j]
        return p

    def accumulate(probabilities):
        total=0
        cumulative=[]
        for p in probabilities:
            total+=p
            cumulative.append(total)
        return cumulative

    def generate(probabilities):
        accumulated=accumulate(probabilities)
        rr=accumulated[len(accumulated)-1]*random()
        i=0
        while accumulated[i]<=rr:
            i+=1
        return i

    motifs=[]

    for i in range(t):
        motifs.append(random_kmer(dna[i]))
    bestMotifs=motifs

    trace=[]
    best_score=float_info.max
    for j in range(n):
        i=randint(0,t-1)
        profile=Profile(dropOneMotif(motifs,i))
        motif_index=generate([probability(dna[i][ll:ll+k],profile)\
                              for ll in range(len(dna[i])-k)])
        motifs[i]=dna[i][motif_index:motif_index+k]
        sc=score(k,motifs)
        if  sc< best_score:
            best_score=sc
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

        def test_ba2a(self):
            '''BA2A 	Implement MotifEnumeration'''
            self.assertEqual('ATA ATT GTT TTT',
                             enumerateMotifs(3,
                                             1,
                                             ['ATTTGGC',
                                            'TGCCTTA',
                                            'CGGTATC',
                                            'GAAAATT']))

        def test_ba2b(self):
            '''BA2B 	Find a Median String'''
            self.assertEqual('GAC',
                             medianString(3,[
                                 'AAATTGACGCAT',
                                 'GACGACCACGTT',
                                 'CGTCAGCGCCTG',
                                 'GCTGAGCACCGG',
                                 'AGTACGGGACAG'
                             ]))


        def test_ba2c(self):
            '''BA2C 	Find a Profile-most Probable k-mer in a String'''
            self.assertEqual(
                'CCGAG',
                mostProbable(
                    'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT',\
                    5,
                    [[0.2, 0.2, 0.3, 0.2, 0.3],
                     [0.4, 0.3, 0.1, 0.5, 0.1],
                     [0.3, 0.3, 0.5, 0.2, 0.4],
                     [0.1, 0.2, 0.1, 0.1, 0.2]]))

        def test_ba2d(self):
            ''' BA2D 	Implement GreedyMotifSearch'''
            motifs=greedyMotifSearch(3,
                                     5,
                                     [
                                         'GGCGTTCAGGCA',
                                         'AAGAATCAGTCA',
                                         'CAAGGAGTTCGC',
                                         'CACGTCAATCAC',
                                         'CAATAATATTCG'
                                     ])
            self.assertEqual(['CAA','CAA','CAA','CAG','CAG'],sorted(motifs))


        def test_ba2e(self):
            '''BA2E 	Implement GreedyMotifSearch with Pseudocounts'''
            motifs=greedyMotifSearch(3,
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
            '''BA2G 	Implement GibbsSampler'''
            x=gibbs(8, 5, 100,[
                'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'])


        def test_ba2h(self):
            '''BA2H 	Implement DistanceBetweenPatternAndStrings'''
            self.assertEqual(
                5,
                distanceBetweenPatternAndStrings('AAA',
                                                 ['TTACCTTAAC',
                                                  'GATATCTGTC',
                                                  'ACGGCGTTCG',
                                                  'CCCTAAAGAG',
                                                  'CGTCAGAGGT']))
        pass



    main()
