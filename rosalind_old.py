# Copyright (C) 2015-2019 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# This file contains a collection of functions to solve the problems
# at rosalind.info.

import helpers as rh, reference_tables as rrt, fasta as f, math, sys, numpy,random, functools

class RosalindException(Exception):
    def __init__( self, message ):
        Exception.__init__(self, message)
        
## Rosalind functions

### How do we sequence antibodies?



# BA4B	Find Substrings of a Genome Encoding a Given Amino Acid String
#
# There are three different ways to divide a DNA string into codons for
# translation, one starting at each of the first three starting positions of
# the string. These different ways of dividing a DNA string into codons are
# called reading frames. Since DNA is double-stranded, a genome has six reading
# frames (three on each strand).
#
# We say that a DNA string Pattern encodes an amino acid string Peptide if
# the RNA string transcribed from either Pattern or its reverse complement
# Pattern translates into Peptide.
#
# Input: A DNA string Text and an amino acid string Peptide.
#
# Return: All substrings of Text encoding Peptide (if any such substrings exist)

def findEncodings(text,peptide):
    def encodes(dna):
        try:
            return prot(dna_to_rna(dna))==peptide
        except KeyError:
            return False    
    candidates=[text[i:i+3*len(peptide)] for i in range(len(text)-3)]
    return [rna for rna in candidates if encodes(rna) or encodes(revc(rna))]
            
# BA4C	Generate the Theoretical Spectrum of a Cyclic Peptide
#
# The workhorse of peptide sequencing is the mass spectrometer, an expensive
# molecular scale that shatters molecules into pieces and then weighs the 
# resulting fragments. The mass spectrometer measures the mass of a molecule
# in daltons (Da); 1 Da is approximately equal to the mass of a single nuclear
# particle (i.e., a proton or neutron).
#
# We will approximate the mass of a molecule by simply adding the number of
# protons and neutrons found in the molecule’s constituent atoms, which yields
# the molecule’s integer mass. For example, the amino acid "Gly", which has 
# chemical formula C2H3ON, has an integer mass of 57,
# since 2·12 + 3·1 + 1·16 + 1·14 = 57. Yet 1 Da is not exactly equal to the mass
# of a proton/neutron, and we may need to account for different naturally
# occurring isotopes of each atom when weighing a molecule. As a result, 
# amino acids typically have non-integer masses (e.g., "Gly" has total mass
# equal to approximately 57.02 Da); for simplicity, however, we will work with
# the integer mass table given in Figure 1.
#
# The theoretical spectrum of a cyclic peptide Peptide, denoted
# Cyclospectrum(Peptide), is the collection of all of the masses of its
# subpeptides, in addition to the mass 0 and the mass of the entire peptide.
# We will assume that the theoretical spectrum can contain duplicate elements,
# as is the case for "NQEL", where "NQ" and "EL" have the same mass.
#
# Input: An amino acid string Peptide.
#
# Return: Cyclospectrum(Peptide).

def cycloSpectrum(peptide,mass=rrt.integer_masses):
    def get_pairs(index_range):
        n=len(index_range)
        return [(i,j) for i in index_range for j in range(i,i+n) if j!=i]
    augmented_peptide=peptide+peptide   # allows easy extraction of substrings
                                        # fromcyclic peptide
    spectrum=[rh.get_mass(augmented_peptide[a:b],mass)\
              for (a,b) in get_pairs(range(len(peptide)))]
    spectrum.append(rh.get_mass('',mass))
    spectrum.append(rh.get_mass(peptide,mass)) 
    spectrum.sort()
    return spectrum

# BA4D	Compute the Number of Peptides of Given Total Mass
#
# In Generate the Theoretical Spectrum of a Cyclic Peptide, we generated the
# theoretical spectrum of a known cyclic peptide. Although this task is
# relatively easy, our aim in mass spectrometry is to solve the reverse problem:
# we must reconstruct an unknown peptide from its experimental spectrum. 
# We will start by assuming that a biologist is lucky enough to generate an
# ideal experimental spectrum Spectrum, which is one coinciding with the
# peptide’s theoretical spectrum. Can we reconstruct a peptide whose
# theoretical spectrum is Spectrum?
#
# Denote the total mass of an amino acid string Peptide as Mass(Peptide).
# In mass spectrometry experiments, whereas the peptide that generated a 
# spectrum is unknown, the peptide’s mass is typically known and is denoted
# ParentMass(Spectrum). Of course, given an ideal experimental spectrum,
# Mass(Peptide) is given by the largest mass in the spectrum.
#
# A brute force approach to reconstructing a peptide from its theoretical
# spectrum would generate all possible peptides whose mass is equal to 
# ParentMass(Spectrum) and then check which of these peptides has theoretical
# spectra matching Spectrum. However, we should be concerned about the running
# time of such an approach: how many peptides are there having mass equal
# to ParentMass(Spectrum)?
#
# Input: An integer m.
#
# Return: The number of linear peptides having integer mass m.
#
# NB, treat peptide as a vector of masses, so amino acids with the same
# mass are the same

def count_peptides_linear(total_mass):
    cache=[]
    masses=list(set(rrt.integer_masses.values()))
    for target_mass in range(total_mass+1):
        total=0
        for amino_acid_mass in masses:
            residual_mass=target_mass-amino_acid_mass
            if residual_mass==0:
                total+=1
            elif residual_mass>0:
                total+=cache[residual_mass]
        cache.append(total)
        
    return total


def get_weight(peptide):
    return sum(rrt.amino_acids[amino_acid].mon_mass for amino_acid in peptide)

def expand(peptides,masses):
    return [peptide+[mass] for peptide in peptides for mass in masses]
 
def mass(peptide):
    return sum([weight for weight in peptide])

def parentMass(spectrum):
    return max(spectrum)

# BA4E 	Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum
#
# Input: A collection of (possibly repeated) integers Spectrum corresponding
#       to an ideal experimental spectrum.
#
# Return: An amino acid string Peptide such that Cyclospectrum(Peptide) = 
#        Spectrum (if such a string exists).

def find_cyclopeptide_sequence(spectrum):
        
    def consistent(peptide,spectrum):
        def count(element,spect):
            return len ([s for s in spect if s==element])
        peptide_spectrum=linearSpectrum(peptide)
        for element in peptide_spectrum:
            if count(element,peptide_spectrum)>count(element,spectrum):
                return False
        return True
    
    def linearSpectrum(peptide):
        def get_pairs():
            return [(i,j) for i in range(len(peptide)) for j in range(i+1,len(peptide))]
        
        result=[sum(peptide[a:b]) for (a,b) in get_pairs()]
        result.append(0)
        result.append(sum(peptide)) 
        result.sort()
        return result 
        
    def cycloSpectrum(peptide):
        def get_pairs(index_range):
            n=len(index_range)
            return [(i,j) for i in index_range for j in range(i,i+n) if j!=i]
        augmented_peptide=peptide+peptide
        result=[sum(augmented_peptide[a:b]) for (a,b) in get_pairs(range(len(peptide)))]
        result.append(0)
        result.append(sum(peptide)) 
        result.sort()
        return result  
      
    peptides=[[]]
    output=[]
    masses=list(set(rrt.integer_masses.values()))
    
    while len(peptides)>0:
        next_peptides=[]
        for peptide in expand(peptides,masses):
            if mass(peptide) == parentMass(spectrum):
                if cycloSpectrum(peptide) == spectrum:
                    output.append(peptide)
            else:
                if consistent(peptide,spectrum):
                    next_peptides.append(peptide)    
        peptides=next_peptides
    return output

# BA4F 	Compute the Score of a Cyclic Peptide Against a Spectrum 
def score(peptide,spectrum,spect_from_peptide=cycloSpectrum):
    return rh.countMatchesInSpectra(spect_from_peptide(peptide),spectrum)

# BA4G 	Implement LeaderboardCyclopeptideSequencing 
def leaderPeptide(n,                                            \
                  spectrum,                                     \
                  masses=list(set(rrt.integer_masses.values())),\
                  spect1=rh.linearSpectrum,                     \
                  spect2=rh.linearSpectrum):
   
    
    leaderPeptide=[]
    leaderBoard=[leaderPeptide]
    
    while len(leaderBoard)>0:
        newBoard=[]
        for peptide in expand(leaderBoard,masses):
            if mass(peptide)==parentMass(spectrum):
                if score(peptide,spectrum,spect2)>\
                   score(leaderPeptide,spectrum,spect2):
                    leaderPeptide=peptide
                newBoard.append(peptide)
            elif mass(peptide)>parentMass(spectrum):
                pass #peptide will be dropped from leader board
            else:
                newBoard.append(peptide)
        leaderBoard=trim(newBoard, spectrum, n,spect1)
    return leaderPeptide

# BA4H 	Generate the Convolution of a Spectrum 
def convolution (spectrum):
    def create_counts(diffs):
        counts={}
        for diff in diffs:
            if not diff in counts:
                counts[diff]=0
            counts[diff]+=1
        return counts        
    diffs=[abs(spectrum[i]-spectrum[j])  \
           for i in range(len(spectrum))  \
           for j in range(i+1,len(spectrum)) if spectrum[i]!=spectrum[j]]
    return sorted([(diff,count)                                      \
                   for diff,count in create_counts(diffs).items()],  \
                  key=lambda x: (-x[1], x[0]))

# BA4I 	Implement ConvolutionCyclopeptideSequencing
#
# Given: An integer M, an integer N, and a collection of
# (possibly repeated) integers Spectrum.
#
# Return: A cyclic peptide LeaderPeptide with amino acids taken only from the
# top M elements (and ties) of the convolution of Spectrum that fall between
# 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).
#
# NB: I had to sort spectrum to pass the testcase in the textbook.

def convolutionCyclopeptideSequencing(m,n,spectrum,low_mass=57,high_mass=200):
    def get_masses_from_spectrum():
        masses=[]
        last_count=0
        for mass,count in convolution(spectrum):
            if low_mass<=mass and mass<=high_mass:
                if len(masses)<m:
                    masses.append(mass)
                    last_count=count
                else:
                    if count==last_count:
                        masses.append(mass)
                    else:
                        return masses
        return masses
    
    return leaderPeptide(n,spectrum,get_masses_from_spectrum(),spect2=rh.cycloSpectrum1)

# BA4L 	Trim a Peptide Leaderboard	
#
# Input: A leaderboard of linear peptides Leaderboard, a linear spectrum
# Spectrum, and an integer N.
#
# Return: The top N peptides from Leaderboard scored against Spectrum.
# Remember to use LinearScore.

def trim(leaderBoard, spectrum,n,spectrum_generator=rh.linearSpectrum):
    if len(leaderBoard)<n:
        return leaderBoard
    peptides_with_scores=[\
        (score(peptide,spectrum,spectrum_generator),peptide)\
        for peptide in leaderBoard]
    peptides_with_scores.sort(reverse=True)   
    (cutoff,_)= peptides_with_scores[n-1]
    return [peptide                                     \
            for (score,peptide) in peptides_with_scores \
            if score>=cutoff]

# Adapter for 'trim', so it will work with peptides as strings

def trim_for_strings(leaderBoard, spectrum,n,\
                     spectrum_generator=rh.linearSpectrum,\
                     masses=rrt.integer_masses):
    numeric_peptides=[]
    trans={}
    for peptide in leaderBoard:
        num=[masses[amino_acid] for amino_acid in peptide]
        numeric_peptides.append(num)
        trans[tuple(num)]=peptide
    trimmed=trim(numeric_peptides, spectrum,n,spectrum_generator)
    return [trans [tuple(peptide)] for peptide in trimmed]

# BA4M 	Solve the Turnpike Problem 
def turnpike(differences):
    def diffs(seq):
        return sorted([x-y for x in seq for y in seq if x>y]) 
    def match(rs,ss):
        for r,s in zip(rs,ss):
            if r!=s:
                return False
        return True    
    def get_length():
        n_zeroes=0
        i=0
        j=len(differences)-1
        index_first_pos=-1
        while i<j:
            if differences[i]+differences[j]!=0:
                raise RosalindException('Mismatch %(i)d and %(j)d'%locals())
            if differences[i]==0:
                n_zeroes+=1
            else:
                index_first_pos=j
            i+=1
            j-=1
        n_zeroes = 2*n_zeroes if len(differences)%2==0 else 2*n_zeroes+1
        if n_zeroes*n_zeroes==len(differences):
            return (index_first_pos,n_zeroes)
        else:
            raise RosalindException('Mismatched lengths')
        

    index_first_pos,length_result =get_length()
    largest=differences[-1]
    indices=[]
    while len(indices)<length_result-2:
        indices.append(index_first_pos+1)
    
    while True:
        result=[0]
        for index in indices:
            result.append(differences[index])
        result.append(largest)
       
        if match(diffs(result),differences[index_first_pos:]):
            return result
        j=len(indices)-1
        while j>0:
            if indices[j]<len(differences)-1:
                indices[j]+=1
                break
            else:
                indices[j]=index_first_pos+1
                j-=1
            
        
# BA4J 	Generate the Theoretical Spectrum of a Linear Peptide 
def linearSpectrumFromString(peptide):
    return rh.linearSpectrum([rrt.integer_masses[a] for a in peptide])

# BA4K 	Compute the Score of a Linear Peptide 
def linearScore(peptide,spectrum):
    return rh.countMatchesInSpectra(linearSpectrumFromString(peptide),spectrum)

def convolution_expanded(spectrum):
    return [diff for (diff,count) in convolution (spectrum) for i in range(count) ]
	  	  	 

### Test cases ###

if __name__=='__main__':
 
    import unittest,fragile
    
    class TestRosalind(unittest.TestCase):
              
 
 

        
        def test_ba6e(self):
            pairs=fragile.find_shared_kmers(3,'AAACTCATC','TTTCAAATC')
            self.assertIn((0, 4),pairs)
            self.assertIn((0, 0),pairs)
            self.assertIn((4, 2),pairs)
            self.assertIn((6, 6),pairs)

            


# BA4B	Find Substrings of a Genome Encoding a Given Amino Acid String
        def test_ba4b(self):
            encodings=findEncodings('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA','MA')
            self.assertEqual(3,len(encodings))
            self.assertIn('ATGGCC',encodings)
            self.assertIn('GGCCAT',encodings)
            self.assertIn('ATGGCC',encodings)

# BA4C	Generate the Theoretical Spectrum of a Cyclic Peptide
        def test_ba4c(self):
            self.assertEqual([0,113,114,128,129,227,242,242,257,355,356,370,371,484],
                             cycloSpectrum('LEQN'))
        
        def test_ba4d(self):
            self.assertEqual(14712706211,count_peptides_linear(1024))
        
        def test_chain(self):
            self.assertAlmostEqual(821.392,get_weight('SKADYEK'),places=3)
        

                
        def test_ba4e(self):
            seq=find_cyclopeptide_sequence([0,113,128,186,241,299,314,427])
            self.assertEqual(6,len(seq))
            self.assertIn([186,128,113],seq)
            self.assertIn([186,113,128],seq)
            self.assertIn([128,186,113],seq)
            self.assertIn([128,113,186],seq)
            self.assertIn([113,186,128],seq)
            self.assertIn([113,128,186],seq)
 
# BA4F 	Compute the Score of a Cyclic Peptide Against a Spectrum 
        def test_ba4f(self):
            self.assertEqual(
                11,
                score(
                    'NQEL',
                    [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]))

# BA4G 	Implement LeaderboardCyclopeptideSequencing 
        def test_ba4g(self):
            self.assertEqual([129, 71, 147, 113],
                             leaderPeptide(
                                 10,
                                 [0, 71, 113, 129, 147, 200, 218,\
                                  260, 313, 331, 347, 389, 460]))            

# BA4H 	Generate the Convolution of a Spectrum 
        def test_ba4h(self):
            self.assertEqual(
                [137, 137, 186, 186, 49, 323],
                convolution_expanded([0, 137, 186, 323]))
            
        #def test_random(self):
            #self.assertEqual([-5.737, -5.217, -5.263, -5.360, -5.958, -6.628, -7.009],
                             #random_genome('CGATACAA',[0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783]))


# BA4I 	Implement ConvolutionCyclopeptideSequencing 
        def test_ba4i(self):
            self.assertEqual([99,71,137,57,72,57],
                             convolutionCyclopeptideSequencing(
                                 20,
                                 60,
                                 [57, 57, 71, 99, 129, 137,\
                                  170, 186, 194, 208, 228,\
                                  265, 285, 299, 307, 323,\
                                  356, 364, 394, 422, 493])) 
            
# BA4J 	Generate the Theoretical Spectrum of a Linear Peptide 
        def test_ba4j(self):
            self.assertEqual(
                [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484],
                linearSpectrumFromString('NQEL'))     

# BA4L 	Trim a Peptide Leaderboard
        def test_ba4l(self):
            self.assertEqual(['LAST', 'ALST'],
                             trim_for_strings(['LAST', 'ALST', 'TLLT', 'TQAS'],
                                  [0,71,87,101,113,158,184,188,259,271,372],
                                  2))

# BA4M 	Solve the Turnpike Problem

        def test_ba4m(self):        
            self.assertEqual([0, 2, 4,7, 10],
                             turnpike([
                                 -10,-8, -7, -6, -5,
                                 -4, -3, -3, -2, -2,
                                 0, 0, 0, 0, 0, 
                                 2, 2, 3, 3, 4,
                                 5, 6, 7, 8, 10
                             ]))


# BA5G 	Compute the Edit Distance Between Two Strings 	 	
                         
# BA5H 	Find a Highest-Scoring Fitting Alignment of Two Strings 	 	
                         
# BA5I 	Find a Highest-Scoring Overlap Alignment of Two Strings 	 	
                         
# BA5J 	Align Two Strings Using Affine Gap Penalties 	 	
                         
# BA5K 	Find a Middle Edge in an Alignment Graph in Linear Space 	 	
                         
# BA5L 	Align Two Strings Using Linear Space 	 	
                         
# BA5M 	Find a Highest-Scoring Multiple Sequence Alignment 	 	
                         
# BA5N 	Find a Topological Ordering of a DAG 
        def test_ba5n(self):
            self.assertEqual([5, 4, 1, 2, 3],
                             topological_order({
                                 1 : [2],
                                 2 : [3],
                                 4 : [2],
                                 5 : [3]
                             }))
            
#LGIS 	Longest Increasing Subsequence
        def test_longestIncreasingSubsequence(self):
            (a,d)=longestIncreasingSubsequence(5,[5, 1, 4, 2, 3])
            self.assertEqual([1,2,3],a)
            self.assertEqual([5,4,3],d)
            (a,d)=longestIncreasingSubsequence(9,[8, 2, 1, 6, 5, 7, 4, 3, 9])
            self.assertEqual([1, 5, 7, 9],a)
            self.assertEqual([8, 6, 5, 4, 3],d)
        
        def test_orf(self):
            string='''>Rosalind_99
            AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'''
            fasta=f.FastaContent(string.split('\n'))            
            peptides=get_reading_frames(fasta)
            self.assertIn('MLLGSFRLIPKETLIQVAGSSPCNLS',peptides)
            self.assertIn('M',peptides)
            self.assertIn('MGMTPRLGLESLLE',peptides)
            self.assertIn('MTPRLGLESLLE',peptides)            
            self.assertEqual(4,len(peptides))
            #for peptide in get_reading_frames(f.FastaFile('./data/rosalind_orf(3).txt') ):
                #print (peptide)

        #def test_wfmd(self):
            #self.assertAlmostEqual(0.772, wfmd(4,6,2,1),3)
            
        #def test_q5(self):
            #for peptide in ['MTAI',
                            #'MLAT',
                            #'MAIT',
                            #'MTAL',
                            #'TALM',
                            #'TLAM']:
                #sp=cycloSpectrum(peptide)
                #if sp==[0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]:
                    #print (peptide, sp)
         
        #def test_q6(self):
            ## 0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333
            #spectrum=[0 ,71, 99, 101, 103, 128, 129, 199, 200, 204, 227, 230, 231,\
                #298, 303, 328, 330, 332, 333]
            #print (spectrum)
            #for peptide in ['TVQ',
                            #'CTV',
                            #'AVQ',
                            #'AQV',
                            #'TCE',
                            #'VAQ']:
                #masses=[rrt.integer_masses[a] for a in peptide]
                #masses.sort()
                #ls=rh.linearSpectrum(masses)
                #if rh.consistent(masses,spectrum):
                    #print (peptide)
                #else:
                    #print (peptide,masses,ls)
                    
    unittest.main(exit=False)
