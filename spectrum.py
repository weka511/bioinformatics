#!/usr/bin/env python

#    Copyright (C) 2019-2024 Simon Crase
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

'''
    Code for Chapters 4 and 11, and utilities for mass spectroscopy
'''

from sys import float_info
from unittest import TestCase, main, skip
import numpy as np

from reference_tables import amino_acids, integer_masses, codon_table, test_masses
from bisect import bisect
from rosalind import dna_to_rna, revc, triplets
from fasta import FastaContent

def get_mass(peptide,masses=integer_masses):
    '''Total mass of amino acids in peptide'''
    return sum([masses[amino_acid] for amino_acid in peptide])

def invert(masses=integer_masses):
    '''
    invert

    Input: A dictionary containing masses of amino acids

    Return: A dictionary mapping masses to amino acids
    '''
    return {v:k for k,v in masses.items()}

def SpectrumGraph(spectrum):
    '''
    BA11A Construct the Graph of a Spectrum

    We represent the masses in a spectrum as a sequence Spectrum of integers  s1,…,sm in increasing order,
    where s1 is zero and sm is the total mass of the (unknown) peptide. We define a labeled graph
    Graph(Spectrum) by forming a node for each element of Spectrum, then connecting nodes si and sj
    by a directed edge labeled by an amino acid a if sj−si is equal to the mass of a.

    As we assumed when sequencing antibiotics, we do not distinguish between amino acids having the same
    integer masses (i.e., the pairs K/Q and I/L).
    '''
    def add_one_point(product,inverted,index=-1):
        '''Add one point to graph'''
        value = spectrum[index] if index>-1 else 0
        for j in range(index+1,len(spectrum)):
            diff = spectrum[j] - value
            if diff in inverted:
                if not value in product:
                    product[value] = []
                for protein in inverted[diff]:
                    product[value].append((spectrum[j],protein))

    inverted = invert(integer_masses)
    product = {}
    add_one_point(product,inverted)
    for i in range(len(spectrum)):
        add_one_point(product,inverted,index=i)
    return product


def DecodeIdealSpectrum(Spectrum):
    '''
    DecodeIdealSpectrum

    Reconstruct a peptide from its ideal spectrum.

     Input: A  list of integers, Spectrum.

     Return: An amino acid string with an ideal spectrum that matches Spectrum.
    '''
    def bfs(adj,paths = [[(0,'')]]):
        closed_paths = []
        while True:
            new_paths=[]
            for path in paths:
                a,_ = path[-1]
                if a in adj:
                    for b,c in adj[a]:
                        new_paths.append(path + [(b,c)])
                else:
                    closed_paths.append(path)
            if len(new_paths)==0:
                return paths + closed_paths
            else:
                paths = new_paths

    for path in bfs(SpectrumGraph(Spectrum)):
        peptide = ''.join(s for (_,s) in path)
        generated_spectrum = linearSpectrum(peptide)
        if generated_spectrum[1:]==Spectrum:
            return peptide

    return None


def create_extended():
    '''
    create_extended
    Extend list of amino acids for testing
    '''
    extended_masses = {'X':4,'Z':5}
    extended_masses.update(integer_masses)
    return extended_masses

def CreatePeptideVector(peptide):
    '''
    CreatePeptideVector

     Convert a Peptide into a Peptide Vector

     Convert a peptide into a binary peptide vector.

     Given an amino acid string Peptide = a1 . . . an of length n, we will represent its
     prefix masses using a binary peptide vector Peptide' with mass(Peptide) coordinates.
     This vector contains a 1 at each of the n prefix coordinates

     mass(a1), mass(a1 a2), . . . , mass(a1 a2 . . . an ) ,
     and it contains a 0 in each of the remaining noise coordinates.

     Input: A peptide P.

     Return: The peptide vector of P.

     Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass
     table for the regular twenty amino acids. Examples sometimes use imaginary amino
     acids X and Z having respective integer masses 4 and 5.
    '''
    extended_masses = create_extended()
    masses = [extended_masses[p] for p in peptide]
    result = []
    for m in masses:
        result = result + ([0]*(m-1))
        result.append(1)
    return result


def CreatePeptide(vector):
    '''
    CreatePeptide

    Convert a Peptide Vector into a Peptide

    '''
    extended_masses = create_extended()
    masses_offset   = [i+1 for i in range(len(vector)) if vector[i]>0]
    masses          = [b-a for (a,b) in zip([0]+masses_offset[:-1],masses_offset)]
    inverted_masses = invert(extended_masses)
    return ''.join( [str(inverted_masses[m][0]) for m in masses])



def conv(S,T,eps=0.001):
    '''
    conv

     Comparing Spectra with the Spectral Convolution

     Comparing Spectraclick to collapse

     Suppose you have two mass spectra, and you want to check if they both were obtained from the same protein;
     you will need some notion of spectra similarity. The simplest possible metric would be to count the number
     of peaks in the mass spectrum that the spectra share, called the shared peaks count;
     its analogue for simplified spectra is the number of masses that the two spectra have in common.

     The shared peaks count can be useful in the simplest cases, but it does not help us if, for example,
     one spectrum corresponds to a peptide contained inside of another peptide from which the second
     spectrum was obtained. In this case, the two spectra are very similar, but the shared peaks count
     will be very small. However, if we shift one spectrum to the right or left, then shared peaks will align.
     In the case of simplified spectra, this means that there is some shift value `x` such that adding
     x to the weight of every element in one spectrum should create a large number of matches in the other spectrum.

     Inputs: Two multisets of positive real numbers S1 and S2

     The size of each multiset is at most 200.

     Return: The largest multiplicity of S1- S2, as well as the absolute value of the number x
     maximizing (S1-S2)(x) (you may return any such value if multiple solutions exist).
    '''
    minkowski_diff = sorted([s-t for s in S for t in T])
    accumulated    = []
    latest         = minkowski_diff[0]
    count          = 1

    for term in minkowski_diff[1:]+[666]:
        if abs(term-latest)<eps:
            count+=1
        else:
            accumulated.append((count,latest))
            latest = term
            count  = 1

    return accumulated[np.argmax([i for i,_ in accumulated])]

def create_lookup(amino_acids=amino_acids):
    '''
    create_lookup

    Creates a lookup table for amino acid masses
    '''
    pairs = sorted([(abbrev,value.mon_mass) for abbrev,value in amino_acids.items()],
                   key =lambda x:x[1])
    pairs.append(('?',float_info.max))
    masses = [mass for (_,mass) in pairs]
    return masses,pairs

def get_abbrev(diff,masses,pairs):
    '''
    get_abbrev

    Find amino acid whose mass best matches a specified difference
    '''
    index = bisect(masses,diff)
    m1 = masses[index]
    m0 = masses[(index-1) if index>0 else 0]
    if index>0 and diff-m0 < m1-diff:
        index-=1
    abbrev,_ = pairs[index]
    return abbrev

def spectrum2protein(ms):
    '''
    spectrum2protein

    spec Inferring Protein from Spectrum

    Introduction to Mass Spectrometry

    In "Calculating Protein Mass", we briefly mentioned an analytic chemical method called mass spectrometry,
    which aims to measure the mass-to-charge ratio of a particle or a molecule. In a mass spectrometer,
    a sample is vaporized (turned into gas), and then particles from the sample are ionized. The resulting
    ions are placed into an electromagnetic field, which separates them based on their charge and mass.
    The output of the mass spectrometer is a mass spectrum, or a plot of ions' possible mass-to-charge ratio
    values with the intensity (actual observed frequency) of ions having these mass-to-charge values.

    For the moment, we will ignore charge and consider a list of the ions' monoisotopic masses as a
    simplified spectrum. Researchers do not possess cheap technology to go in and examine a protein one
    amino acid at a time (molecules are too submicroscopic). Instead, to determine a protein's structure,
    we will split several copies of the protein into smaller pieces, then weigh the resulting fragments.
    To do this, we assume that each cut (breakage point) occurs between two amino acids and that we can
    measure the mass of the resulting pieces for all possible cuts.

    For example, the (unknown) protein "PRTEIN" can be cut in five possible ways:
    "P" and "RTEIN"; "PR" and "TEIN"; "PRT" and "EIN"; "PRTE" and "IN"; "PRTEI" and "N".
    We then can measure the masses of all fragments, including the entire string. The "left" end of
    a protein is called its N-terminus, and the ions corresponding to the protein string's prefixes
    (P, PR, PRT, PRTE, PRTEI) are called b-ions. The "right" end of the protein is called its C-terminus,
    and the ions corresponding to the string's suffixes (N, IN, EIN, TEIN, RTEIN) are called y-ions.
    The difference in the masses of two adjacent b-ions (or y-ions) gives the mass of one amino acid residue;
    for example, the difference between the masses of "PRT" and "PR" must be the mass of "T." By extension,
    knowing the masses of every b-ion of a protein allows us to deduce the protein's identity.

    The prefix spectrum of a weighted string is the collection of all its prefix weights.

    Input: A list L of n (n<=100) positive real numbers.

    Return: A protein string of length n-1 whose prefix spectrum is equal to L (if multiple solutions exist,
    you may output any one of them). Consult the monoisotopic mass table.
    '''
    masses,pairs = create_lookup()
    return ''.join([get_abbrev(diff,masses,pairs) for diff in [m1-m0 for m0,m1 in zip(ms[:-1],ms[1:])]])

def complete_spectrum(P):
    def spectrum(S):
        return sum([amino_acids[s].mon_mass for s in S])
    prefixes = [P[:i] for i in range(1,len(P))] +[P]
    suffixes = [P[i:] for i in range(1,len(P))]
    ss= [spectrum(p) for p in prefixes + suffixes]
    return ss

def prsm(s,R):
    '''
    prsm

    Match a Spectrum to a Protein

    Searching the Protein Database

    Many proteins have already been identified for a wide variety of organisms. Accordingly,
    there are a large number of protein databases available, and so the first step after
    creating a mass spectrum for an unidentified protein is to search through these databases
    for a known protein with a highly similar spectrum. In this manner, many similar proteins
    found in different species have been identified, which aids researchers in determining protein function.

    In "Comparing Spectra with the Spectral Convolution", we introduced the spectral convolution
    and used it to measure the similarity of simplified spectra. In this problem, we would like
    to extend this idea to find the most similar protein in a database to a spectrum taken from
    an unknown protein. Our plan is to use the spectral convolution to find the largest possible
    number of masses that each database protein shares with our candidate protein after shifting,
    and then select the database protein having the largest such number of shared masses.

    Inputs: A positive integer n followed by a collection of n protein strings s1, s2, ..., sn and a multiset R
            of positive numbers (corresponding to the complete spectrum of some unknown protein string).

            Return: The maximum multiplicity of R-S[sk]
                    taken over all strings sk, followed by the string sk for which this
                    maximum multiplicity occurs (you may output any such value if multiple solutions exist).
    '''
    def count(c):
        return c[1][0]
    m = 0
    i = 0
    Ss = [(P,complete_spectrum(P)) for P in s]

    Cs = [(P,conv(R,S1,eps=0.00001)) for (P,S1) in Ss]

    index = np.argmax([count(c) for c in Cs])

    _,(b,_) = Cs[index]

    return b,s[index]



def full(L,epsilon=0.000001):
    '''

    Inferring Peptide from Full Spectrum

    Inputs : A list L containing 2n+3 positive real numbers (n<=100).
             The first number in L is the parent mass of a peptide P, and all
             other numbers represent the masses of some b-ions and y-ions of P
             (in no particular order). You may assume that if the mass of a b-ion is present,
             then so is that of its complementary y-ion, and vice-versa.

    Return: A protein string t
            of length n for which there exist two positive real numbers w1 and w2
            such that for every prefix p and suffix s of t, each of w(p)+w1 and w(s)+w2
            is equal to an element of L.
            (In other words, there exists a protein string whose t-prefix and
            t-suffix weights correspond to the non-parent mass values of L.)
            If multiple solutions exist, you may output any one.
    '''

    def get_n():
        '''
        Get n, and verify that it is odd
        '''
        n = (len(L)-3)//2
        assert(2*n+3==len(L))
        return n

    def extract(seq,candidates):
        '''
         Extract prefixes or suffixes
         Inputs: seq  Start of prefixes of suffixes
                 candidates

         Returns: prefixes or suffixes
        '''
        while True:
            key = seq[-1]
            if not key in candidates:
                return seq
            succs = candidates[key]
            _,j,_,_,_,_,_ = succs[0]
            seq.append(j)

    masses,pairs   = create_lookup()
    n              = get_n()

    # Idea: L = prefixes + suffixes, and each of prefexes and suffixes will have many elements whose
    #           differences are equal to the mass of one amino acid. Start by compouting all differences

    diffs          = [(i,j,L[i],L[j],abs(L[i] - L[j])) for i in range(1,len(L)) for j in range(i+1,len(L))]

    # Now compute a collection of candidates, i.e. differences that are close to the mass py one amino acid.

    candidates     = {}
    for i,j,l1,l2,diff in diffs:
        abbrev    = get_abbrev(diff,masses,pairs)
        candidate = abbrev if abs(diff-amino_acids[abbrev].mon_mass)<epsilon else None
        if candidate != None:
            if not i in candidates:
                candidates[i]=[]
            candidates[i].append((i,j,l1,l2,diff,candidate,abs(diff-amino_acids[abbrev].mon_mass)))

    # Now support the data so that they are organized in ascending order by difference from nearest amino acid
    for key in candidates:
        candidates[key]= sorted(candidates[key],key=lambda x:x[6])

    prefixes = extract( [min(candidates.keys())],candidates)
    suffixes = extract([min([r for r in candidates.keys() if not r in prefixes])],candidates)
    # Our algorithm is a bit greedy, so there may be some overlap between prefixes and suffices
    # Thius must be fixed!

    while len(prefixes)>len(suffixes):
        for r in suffixes:
            if r in prefixes:
                ii = prefixes.index(r)
                prefixes = prefixes[:ii] + prefixes[ii+1:]

    return ''.join([candidates[l][0][5] for l in prefixes][:-1])


def sgra(L,Alphabet=amino_acids,epsilon=0.001):
    '''
     sgra

     Using the Spectrum Graph to Infer Peptides

     In this problem, we say that a weighted string s=s1s2...sn
     matches L if there is some increasing sequence of positive real numbers (w1,w2,,,,,wn+1)
     in L such that w(s1)=w2-w1, w(s2)=w3-w2, ..., and w(sn)=wn+1-wn

     Input: A list L (of length at most 100) containing positive real numbers.

     Return: The longest protein string that matches the spectrum graph of L (if multiple solutions exist,
      you may output any one of them). Consult the monoisotopic mass table.

     NB: we want the longest path trough the spectrum graph

    '''

    def create_spectrum_graph():
        ''''
        create_spectrum_graph

        For a weighted alphabet A and a collection L of positive real numbers, the spectrum graph of L
        is a digraph constructed in the following way. First, create a node for every real number in L.
        Then, connect a pair of nodes with a directed edge (u,v) if v>u and v-u is equal to the weight of a single symbol in A

        We may then label the edge with this symbol.
        '''
        masses,pairs = create_lookup()
        G            = {}
        for u in L:
            for v in L:
                if u<v:
                    abbrev = get_abbrev(v-u,masses,pairs)
                    if abs(v-u-amino_acids[abbrev].mon_mass)<epsilon:
                        if not u in G:
                            G[u]=[]
                        G[u].append((abbrev,v))
        return G

    def dfs(key,G,path):
        '''
        dfs

        Depth first search to build up peptides

        '''
        Runs.append(path)
        if key in Outs: return
        if not key in G: return
        for amino_acid,mass in G[key]:
            dfs(mass,G,path+[amino_acid])
            Outs.add(mass)

    G    = create_spectrum_graph()
    Outs = set()     # Used by dfs to determine whter it has already processed node
    Runs = []        # The peptids string that are being built by dfs

    for key in sorted(G.keys()):
        dfs(key,G,[])

    return ''.join(Runs[np.argmax([len(run) for run in Runs])])


def Turnpike(D,check=False):
    '''
       Turnpike

    BA4M 	Solve the Turnpike Problem

    Solve the Turnpike Problem
    Parameters:
        D     The differences
        check Indicates whether differences are to be checked

    Based Mark Weiss's treatment
        https://users.cs.fiu.edu/~weiss/cop3337_f99/assignments/turnpike.pdf

    '''

    def find_remaining_points(X,D,first,last):
        '''
        Extend a partial solution by adding more points.
          Parameters
              X         Partial solution - typically contains NaNs for unassigned points,
                                           with known values filled in from the ends
              D
              first     Index of end of knownn values at beginning X
              last      Index of start of know value at end of X
        '''
        def get_set_diffs(diffs):
            '''
            Ascertain whether a some partial set of differences is a subset of D
            '''
            diffs.sort()
            set_diffs=[]
            i=0
            for diff in diffs:
                found=False
                while i<len(D):
                    if  D[i]<diff:
                        set_diffs.append(D[i])
                    elif D[i]==diff:
                        found=True
                        i+=1
                        break
                    i+=1
                if not found:
                    return None
            while i<len(D):
                set_diffs.append(D[i])
                i+=1
            return set_diffs

        def explore(candidate,X,first,last):
            '''
            Explore lower levels of tree
            '''
            # Constuct set of differences between candidate and known members of X
            diffs=[abs(candidate-x) for x in X if not np.isnan(x) and x!=candidate]
            set_diffs=get_set_diffs(diffs)
            if set_diffs==None:
                return None
            elif len(set_diffs)==0:
                return X
            else:
                return find_remaining_points(X,set_diffs,first,last)

        # There are two cases to consider: either the largest remaining unprocessed value in D
        # is part of the solution, or it isn't. We will explore these two cases separately.
        # We maintain a tree of data structures, so we can explore the tree of solutions
        # see https://users.cs.fiu.edu/~weiss/cop3337_f99/assignments/turnpike.pdf for details
        x_max          = D[-1]
        XX             = X[:]           # Clone this so the
        XX[last-1]     = x_max  # Add candidate at end
        trial_solution = explore(x_max,XX,first,last-1) #process level below - one fewer unknown at end
        if trial_solution==None:  # largest remaining unprocessed value was a false lead
            XX=X[:]
            XX[first+1]=X[-1]-x_max # Added new candiate at beginning
            return explore(X[-1]-x_max,XX,first+1,last) #process level below - one fewer unknown at start
        else:
            return trial_solution

    def check_diffs(reconstruction):
        '''
        Verify that a particular reconstruction does indeed give rise to original differences
        '''
        diffs=[a-b for a in reconstruction for b in reconstruction]
        diffs.sort()
        if len(diffs)!=len(D):
            raise RosalindException ('Length of reconstructed diffs ({0}) does not match length of D ({1}) '.format(len(diffs), len(D)))
        mismatches=0
        for a,b in zip(D,diffs):
            if a!=b:
                mismatches+=1
                print (a,b)

        if mismatches>0:
            raise RosalindException ('Found {0} mismatches'.format(mismatches))
        return diffs

    # Start by initializing array of points. We know that its length
    # must be the square root of the array of differences.
    len_D=len (D)
    len_X=int(np.sqrt(len_D))
    X=[float('nan')]*len_X    #We fill in all values as "unknown"
    X[0]=0                    #Actually we are given the first value, zero
    X[-1]=D[-1]               # We also know that the last point must match the last difference.
    reconstruction= find_remaining_points(X,[d for d in D[:-1] if d>0],0,-1)
    if check:
        check_diffs(reconstruction)
    return reconstruction




def prot(rna,table=codon_table):
    '''
     BA4A	Translate an RNA String into an Amino Acid String
     PROT Translating RNA into Protein

     The 20 commonly occurring amino acids are abbreviated by using 20 letters from
     the Roman alphabet (all letters except for B, J, O, U, X, and Z). Protein strings
     are constructed from these 20 symbols. Henceforth, the term genetic string will
     incorporate protein strings along with DNA strings and RNA strings.

     The RNA codon table dictates the details regarding the encoding of specific
     codons into the amino acid alphabet.

     Input: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).

     Return: The protein string encoded by s.

     NB: I have allowed an extra parameter to deal with alternatives, such as the
     Mitochondrial codes (http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
     '''
    return ''.join([table[codon] for codon in triplets(rna) if table[codon]!=';'])

def findEncodings(text,peptide):
    '''
    BA4B	Find Substrings of a Genome Encoding a Given Amino Acid String

     There are three different ways to divide a DNA string into codons for
     translation, one starting at each of the first three starting positions of
     the string. These different ways of dividing a DNA string into codons are
     called reading frames. Since DNA is double-stranded, a genome has six reading
     frames (three on each strand).

     We say that a DNA string Pattern encodes an amino acid string Peptide if
     the RNA string transcribed from either Pattern or its reverse complement
     Pattern translates into Peptide.

     Input: A DNA string Text and an amino acid string Peptide.

     Return: All substrings of Text encoding Peptide (if any such substrings exist)
    '''
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
    # Input: An amino acid string Peptide (1 char abbreviations).
    #
    # Return: Cyclospectrum(Peptide).

def cycloSpectrum(peptide,mass=integer_masses):

    # get_pairs
    #
    # Inputs: index_range
    #
    # Return:  Pairs of indices delimiting sublists of peptide

    def get_pairs(index_range):
        return [(i,j) for i in index_range for j in range(i,i+len(index_range)) if j!=i]

    augmented_peptide = peptide+peptide   # allows easy extraction of substrings
                                          # fromcyclic peptide
    spectrum          = [get_mass(augmented_peptide[a:b],mass) for (a,b) in get_pairs(range(len(peptide)))]
    spectrum.append(get_mass('',mass))
    spectrum.append(get_mass(peptide,mass))
    spectrum.sort()
    return spectrum

def cycloSpectrum1(peptide):
    def get_pairs(index_range):
        n=len(index_range)
        return [(i,j) for i in index_range for j in range(i,i+n) if j!=i]
    augmented_peptide=peptide+peptide
    result=[sum(augmented_peptide[a:b]) for (a,b) in get_pairs(range(len(peptide)))]
    result.append(0)
    result.append(sum(peptide))
    result.sort()
    return result

def count_peptides_linear(total_mass):
    '''
    BA4D	Compute the Number of Peptides of Given Total Mass

     In Generate the Theoretical Spectrum of a Cyclic Peptide, we generated the
     theoretical spectrum of a known cyclic peptide. Although this task is
     relatively easy, our aim in mass spectrometry is to solve the reverse problem:
     we must reconstruct an unknown peptide from its experimental spectrum.
     We will start by assuming that a biologist is lucky enough to generate an
     ideal experimental spectrum Spectrum, which is one coinciding with the
     peptide’s theoretical spectrum. Can we reconstruct a peptide whose
     theoretical spectrum is Spectrum?

     Denote the total mass of an amino acid string Peptide as Mass(Peptide).
     In mass spectrometry experiments, whereas the peptide that generated a
     spectrum is unknown, the peptide’s mass is typically known and is denoted
     ParentMass(Spectrum). Of course, given an ideal experimental spectrum,
     Mass(Peptide) is given by the largest mass in the spectrum.

     A brute force approach to reconstructing a peptide from its theoretical
     spectrum would generate all possible peptides whose mass is equal to
     ParentMass(Spectrum) and then check which of these peptides has theoretical
     spectra matching Spectrum. However, we should be concerned about the running
     time of such an approach: how many peptides are there having mass equal
     to ParentMass(Spectrum)?

     Input: An integer m.

     Return: The number of linear peptides having integer mass m.

     NB, treat peptide as a vector of masses, so amino acids with the same
     mass are the same
    '''
    cache=[]
    masses=list(set(integer_masses.values()))
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
    '''
     Input:  peptide   List of amino acids

     Return:  Monoisotopic mass of peptide
    '''
    return sum(amino_acids[amino_acid].mon_mass for amino_acid in peptide)

def expand(peptides,masses):
    return [peptide+[mass] for peptide in peptides for mass in masses]

def mass(peptide):
    return sum([weight for weight in peptide])

def parentMass(spectrum):
    return max(spectrum)



def find_cyclopeptide_sequence(spectrum):
    '''
    BA4E 	Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum

     Input: A collection of (possibly repeated) integers Spectrum corresponding
           to an ideal experimental spectrum.

     Return: An amino acid string Peptide such that Cyclospectrum(Peptide) =
            Spectrum (if such a string exists).
    '''
    # isConsistent
    #
    # Determine whether peptide is consistent with spectrum
    def isConsistent(peptide):
        # count
        #
        # Determine number of items in spect that match specified element

        def count(element,spect):
            return len ([s for s in spect if s==element])

        peptide_spectrum = linearSpectrum(peptide)

        for element in peptide_spectrum:
            if count(element,peptide_spectrum)>count(element,spectrum):
                return False
        return True

    # cycloSpectrum
    #
    # Compute spectrum for cyclic peptide
    #
    # Inputs:  peptide   Peptide represented as a list of masses
    #
    # Returns: spectrum of peptide

    def cycloSpectrum(peptide):

        # get_pairs
        #
        # Inputs: index_range
        #
        # Return:  Pairs of indices delimiting sublists of peptide

        def get_pairs(index_range):
            return [(i,j) for i in index_range for j in range(i,i+len(index_range)) if j!=i]

        augmented_peptide = peptide+peptide
        result            = [sum(augmented_peptide[a:b]) for (a,b) in get_pairs(range(len(peptide)))]
        result.append(0)
        result.append(sum(peptide))
        result.sort()
        return result

    peptides = [[]]
    output   = []
    masses   = list(set(integer_masses.values()))

    while len(peptides)>0:
        next_peptides=[]
        for peptide in expand(peptides,masses):
            if mass(peptide) == parentMass(spectrum):
                if cycloSpectrum(peptide) == spectrum:
                    output.append(peptide)
            else:
                if isConsistent(peptide):
                    next_peptides.append(peptide)
        peptides=next_peptides

    return output


def countMatchesInSpectra(spect1,spect2):
    i1=0
    i2=0
    count=0
    while i1<len(spect1) and i2<len(spect2):
        diff=spect1[i1]-spect2[i2]
        if diff==0:
            count+=1
            i1+=1
            i2+=1
        elif diff<0:
            i1+=1
        else:
            i2+=1
    return count

def score(peptide,spectrum,spect_from_peptide=cycloSpectrum):
    '''BA4F 	Compute the Score of a Cyclic Peptide Against a Spectrum'''
    return countMatchesInSpectra(spect_from_peptide(peptide),spectrum)

def linearSpectrum(peptide):
    '''BA4J 	Generate the Theoretical Spectrum of a Linear Peptide'''
    def get_pairs():
        return [(i,j) for i in range(len(peptide)) for j in range(len(peptide)+1) if i<j]
    result = [sum(peptide[i:j]) for (i,j) in get_pairs()]
    result.append(0)
    result.sort()
    return result

def leaderPeptide(n,
                  spectrum,
                  masses=list(set(integer_masses.values())),
                  spect1=linearSpectrum,
                  spect2=linearSpectrum):
    '''BA4G 	Implement LeaderboardCyclopeptideSequencing'''
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

def convolution (spectrum):
    '''BA4H 	Generate the Convolution of a Spectrum'''
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

def convolution_expanded(spectrum):
    ''' BA4H 	Generate the Convolution of a Spectrum'''
    return [diff for (diff,count) in convolution (spectrum) for i in range(count) ]


def convolutionCyclopeptideSequencing(m,n,spectrum,low_mass=57,high_mass=200):
    '''
     BA4I 	Implement ConvolutionCyclopeptideSequencing

     Given: An integer M, an integer N, and a collection of
     (possibly repeated) integers Spectrum.

     Return: A cyclic peptide LeaderPeptide with amino acids taken only from the
     top M elements (and ties) of the convolution of Spectrum that fall between
     57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).

     NB: I had to sort spectrum to pass the testcase in the textbook.
    '''
    def get_masses_from_spectrum():
        masses     = []
        last_count = 0
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

    return leaderPeptide(n,spectrum,get_masses_from_spectrum(),spect2=cycloSpectrum1)



def linearSpectrumFromString(peptide):
    '''BA4J 	Generate the Theoretical Spectrum of a Linear Peptide'''
    return linearSpectrum([integer_masses[a] for a in peptide])

def linearScore(peptide,spectrum):
    '''BA4K 	Compute the Score of a Linear Peptide'''
    return countMatchesInSpectra(linearSpectrumFromString(peptide),spectrum)

def trim(leaderBoard, spectrum,n,spectrum_generator=linearSpectrum):
    '''
    BA4L 	Trim a Peptide Leaderboard

    Input: A leaderboard of linear peptides Leaderboard, a linear spectrum
    Spectrum, and an integer N.

    Return: The top N peptides from Leaderboard scored against Spectrum.
    Remember to use LinearScore.
    '''
    if len(leaderBoard)<n:
        return leaderBoard
    peptides_with_scores=[
        (score(peptide,spectrum,spectrum_generator),peptide)
        for peptide in leaderBoard]
    peptides_with_scores.sort(reverse=True)
    (cutoff,_)= peptides_with_scores[n-1]
    return [peptide
            for (score,peptide) in peptides_with_scores
            if score>=cutoff]

def trim_for_strings(leaderBoard, spectrum,n,
                     spectrum_generator=linearSpectrum,
                     masses=integer_masses):
    '''BA4L Adapter for 'trim', so it will work with peptides as strings'''

    numeric_peptides=[]
    trans={}
    for peptide in leaderBoard:
        num=[masses[amino_acid] for amino_acid in peptide]
        numeric_peptides.append(num)
        trans[tuple(num)]=peptide
    trimmed=trim(numeric_peptides, spectrum,n,spectrum_generator)
    return [trans [tuple(peptide)] for peptide in trimmed]



def splc(fasta):
    '''
    SPLC	RNA Splicing

     After identifying the exons and introns of an RNA string, we only need to
     delete the introns and concatenate the exons to form a new string
     ready for translation.

     Input: A DNA string s (of length at most 1 kbp) and a collection of substrings
            of s acting as introns. All strings are given in FASTA format.

     Return: A protein string resulting from transcribing and translating the
             exons of s. (Note: Only one solution will exist for the dataset provided.)
    '''
    (_,dna)=fasta[0]
    for i in range(1,len(fasta)):
        (l,intron)=fasta[i]
        fragments=dna.split(intron)
        if len(fragments)>1:
            dna=''.join(fragments)
    return prot(dna_to_rna(dna))

def SequencePeptide(spectral, protein_masses = integer_masses):
    '''
    BA11E Sequence a Peptide

    Input: A spectral vector S.

    Returns: A peptide with maximum score against S. For masses with more than one amino acid, any choice may be used.
    '''
    def create_inverse_masses():
        product = {}
        for protein,mass in protein_masses.items():
            if not mass in product:
                product[mass] = []
            product[mass].append(protein)
        return product

    def create_DAG(inverse_masses):
        Nodes = [0] + spectral
        Edges = []
        for i in range(len(Nodes)):
            Edges.append([])
            for j in range(i+1,len(Nodes)):
                if j-i in inverse_masses:
                    Edges[-1].append((j,j-i))
        return Nodes,Edges

    inverse_masses = create_inverse_masses()
    Nodes,Edges = create_DAG(inverse_masses)


    return 'XZZXX'

if __name__=='__main__':

    class TestSpectrum(TestCase):

        def test_spec(self):
            self.assertEqual('WMQS',
                             spectrum2protein([3524.8542,3710.9335,3841.974,3970.0326,4057.0646]))

        def test_conv(self):
            m,x=conv([186.07931, 287.12699 ,548.20532 ,580.18077 ,681.22845, 706.27446, 782.27613 ,968.35544, 968.35544],
                 [101.04768, 158.06914 ,202.09536 ,318.09979 ,419.14747, 463.17369])
            self.assertEqual(3,m)
            self.assertAlmostEqual(85.03163,x,places=5)

        def test_prsm(self):
            m,s_max = prsm(['GSDMQS',
                            'VWICN',
                            'IASWMQS',
                            'PVSMGAD'],
                           [445.17838,
                            115.02694,
                            186.07931,
                            314.13789,
                            317.1198,
                            215.09061]
                           )
            self.assertEqual(3,m)
            self.assertEqual('GSDMQS',s_max)

        def test_full(self):
            self.assertEqual('KEKEP',
                             full([
                                1988.21104821,
                                610.391039105,
                                738.485999105,
                                766.492149105,
                                863.544909105,
                                867.528589105,
                                992.587499105,
                                995.623549105,
                                1120.6824591,
                                1124.6661391,
                                1221.7188991,
                                1249.7250491,
                                1377.8200091
                             ]))

        def test_sgra(self):
            '''SGRA Using the Spectrum Graph to Infer Peptides '''
            self.assertEqual('WMSPG',sgra([
                3524.8542,
                3623.5245,
                3710.9335,
                3841.974,
                3929.00603,
                3970.0326,
                4026.05879,
                4057.0646,
                4083.08025
            ]))


        def test_ba4b(self):
            '''BA4B	Find Substrings of a Genome Encoding a Given Amino Acid String'''
            encodings=findEncodings('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA','MA')
            self.assertEqual(3,len(encodings))
            self.assertIn('ATGGCC',encodings)
            self.assertIn('GGCCAT',encodings)
            self.assertIn('ATGGCC',encodings)

        def test_ba4c(self):
            '''BA4C	Generate the Theoretical Spectrum of a Cyclic Peptide'''
            self.assertEqual([0,113,114,128,129,227,242,242,257,355,356,370,371,484],
                             cycloSpectrum('LEQN'))


        def test_ba4d(self):
            '''BA4D	Compute the Number of Peptides of Given Total Mass'''
            self.assertEqual(14712706211,count_peptides_linear(1024))


        def test_chain(self):
            self.assertAlmostEqual(821.392,get_weight('SKADYEK'),places=3)

        def test_ba4e(self):
            '''BA4E 	Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum'''
            seq = find_cyclopeptide_sequence([0,113,128,186,241,299,314,427])
            self.assertEqual(6,len(seq))
            self.assertIn([186,128,113],seq)
            self.assertIn([186,113,128],seq)
            self.assertIn([128,186,113],seq)
            self.assertIn([128,113,186],seq)
            self.assertIn([113,186,128],seq)
            self.assertIn([113,128,186],seq)

        def test_ba4f(self):
            ''' BA4F 	Compute the Score of a Cyclic Peptide Against a Spectrum'''
            self.assertEqual(
                11,
                score(
                    'NQEL',
                    [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]))

        def test_ba4g(self):
            '''BA4G 	Implement LeaderboardCyclopeptideSequencing'''
            self.assertEqual([129, 71, 147, 113],
                             leaderPeptide(
                                 10,
                                 [0, 71, 113, 129, 147, 200, 218,\
                                  260, 313, 331, 347, 389, 460]))

        def test_ba4h(self):
            ''' BA4H 	Generate the Convolution of a Spectrum'''
            self.assertEqual(
                [137, 137, 186, 186, 49, 323],
                convolution_expanded([0, 137, 186, 323]))

        def test_ba4i(self):
            '''BA4I 	Implement ConvolutionCyclopeptideSequencing'''
            self.assertEqual([99,71,137,57,72,57],
                             convolutionCyclopeptideSequencing(
                                 20,
                                 60,
                                 [57, 57, 71, 99, 129, 137,\
                                  170, 186, 194, 208, 228,\
                                  265, 285, 299, 307, 323,\
                                  356, 364, 394, 422, 493]))


        def test_ba4j(self):
            '''BA4J 	Generate the Theoretical Spectrum of a Linear Peptide'''
            self.assertEqual(
                [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484],
                linearSpectrumFromString('NQEL'))


        def test_ba4l(self):
            '''BA4L 	Trim a Peptide Leaderboard'''
            self.assertEqual(['LAST', 'ALST'],
                             trim_for_strings(['LAST', 'ALST', 'TLLT', 'TQAS'],
                                  [0,71,87,101,113,158,184,188,259,271,372],
                                  2))

        def test_ba4m(self):
            self.assertEqual([0, 3, 6, 8, 10],
                              Turnpike(
                                  [-10, -8, -7, -6, -5,
                                   -4, -3 ,-3 ,-2 ,-2,
                                   0, 0, 0, 0, 0,
                                   2, 2, 3, 3, 4,
                                   5, 6, 7, 8, 10],
                                  check=True))
            # Data from rosalind_ba4m_5_output - my submission dated Sept. 20, 2017, 11:38 a.m.
            self.assertEqual([0, 9, 12, 27, 36, 43, 46, 60, 63, 82, 94, 97, 100, 115, 134, 141, 158, 171, 185, 193,
                              200, 211, 230, 235, 243, 260, 276, 294, 305, 314, 319, 335, 351, 355, 368, 370, 380,
                              385, 395, 403, 412, 415, 429, 432, 449, 466, 471, 488, 491, 509, 527, 532, 540, 548,
                              561, 571, 579, 591, 598, 608, 621, 636, 646, 661, 673, 675, 690, 696, 701, 706, 717,
                              721, 738, 739, 748, 766, 782, 787, 796, 804, 807, 816, 823, 836, 855, 873, 876, 891,
                              905, 907, 924, 937, 952, 969, 973, 992, 1008, 1011, 1024, 1029, 1039, 1057, 1062, 1067,
                              1070, 1078, 1079, 1089, 1106, 1111, 1121, 1136, 1144, 1151, 1170, 1176, 1185, 1190,
                              1192, 1199, 1201, 1210, 1228, 1231, 1234, 1239, 1245, 1258, 1265],
                             Turnpike([
                                 -1265,-1258,-1256,-1253,-1249,-1246,-1245,-1239,-1238,-1236,-1234,-1233,
                                 -1231,-1231,-1230,-1229,-1228,-1227,-1225,-1222,-1222,-1222,-1222,-1219,
                                 -1219,-1219,-1218,-1216,-1215,-1212,-1212,-1210,-1209,-1207,-1205,-1204,
                                 -1203,-1202,-1202,-1201,-1201,-1201,-1199,-1199,-1198,-1198,-1198,-1196,
                                 -1195,-1195,-1193,-1192,-1192,-1192,-1191,-1190,-1190,-1189,-1188,-1188,
                                 -1187,-1185,-1185,-1185,-1185,-1183,-1183,-1183,-1182,-1182,-1181,-1180,
                                 -1179,-1178,-1176,-1176,-1176,-1176,-1174,-1174,-1174,-1173,-1172,-1171,
                                 -1171,-1171,-1170,-1168,-1168,-1168,-1167,-1167,-1165,-1165,-1165,-1165,
                                 -1164,-1164,-1164,-1163,-1163,-1163,-1161,-1161,-1158,-1158,-1158,-1158,
                                 -1157,-1156,-1156,-1155,-1154,-1153,-1152,-1151,-1151,-1150,-1150,-1149,
                                 -1149,-1149,-1149,-1148,-1147,-1147,-1146,-1146,-1145,-1145,-1144,-1144,
                                 -1143,-1143,-1142,-1142,-1142,-1141,-1140,-1140,-1139,-1139,-1139,-1139,
                                 -1138,-1137,-1137,-1136,-1136,-1135,-1134,-1134,-1134,-1134,-1133,-1132,
                                 -1132,-1131,-1131,-1131,-1130,-1130,-1130,-1129,-1128,-1128,-1127,-1127,
                                 -1127,-1125,-1124,-1124,-1124,-1124,-1124,-1124,-1122,-1121,-1119,-1119,
                                 -1117,-1117,-1117,-1116,-1116,-1116,-1115,-1113,-1113,-1113,-1112,-1111,
                                 -1111,-1110,-1110,-1110,-1109,-1109,-1108,-1108,-1108,-1107,-1107,-1107,
                                 -1106,-1105,-1105,-1105,-1104,-1104,-1103,-1102,-1102,-1101,-1101,-1100,
                                 -1100,-1100,-1099,-1099,-1098,-1098,-1098,-1097,-1097,-1096,-1095,-1095,
                                 -1094,-1094,-1094,-1094,-1094,-1093,-1093,-1093,-1092,-1091,-1091,-1090,
                                 -1090,-1090,-1089,-1088,-1088,-1088,-1087,-1087,-1087,-1086,-1085,-1085,
                                 -1084,-1084,-1084,-1082,-1081,-1081,-1080,-1080,-1079,-1079,-1079,-1078,
                                 -1078,-1077,-1077,-1076,-1076,-1076,-1076,-1076,-1075,-1075,-1075,-1074,
                                 -1073,-1073,-1073,-1073,-1072,-1070,-1070,-1070,-1070,-1070,-1070,-1069,
                                 -1069,-1069,-1068,-1068,-1067,-1067,-1067,-1066,-1065,-1065,-1065,-1065,
                                 -1063,-1063,-1062,-1062,-1062,-1061,-1061,-1061,-1060,-1060,-1060,-1060,
                                 -1058,-1058,-1058,-1058,-1058,-1058,-1057,-1057,-1057,-1056,-1055,-1055,
                                 -1054,-1054,-1054,-1054,-1053,-1053,-1052,-1052,-1052,-1051,-1051,-1051,
                                 -1051,-1051,-1050,-1050,-1049,-1049,-1048,-1048,-1047,-1047,-1046,-1046,
                                 -1046,-1046,-1045,-1045,-1044,-1044,-1043,-1043,-1043,-1043,-1043,-1043,
                                 -1042,-1042,-1042,-1041,-1041,-1040,-1039,-1039,-1039,-1039,-1039,-1038,
                                 -1036,-1036,-1036,-1036,-1035,-1035,-1035,-1035,-1035,-1034,-1034,-1034,
                                 -1034,-1033,-1032,-1032,-1031,-1031,-1030,-1030,-1030,-1030,-1029,-1029,
                                 -1029,-1029,-1029,-1028,-1028,-1028,-1028,-1027,-1027,-1027,-1027,-1026,
                                 -1026,-1025,-1024,-1024,-1024,-1024,-1024,-1023,-1023,-1022,-1021,-1021,
                                 -1021,-1021,-1021,-1020,-1020,-1019,-1019,-1019,-1018,-1018,-1017,-1017,
                                 -1017,-1017,-1017,-1016,-1016,-1016,-1015,-1015,-1015,-1015,-1014,-1014,
                                 -1014,-1014,-1012,-1012,-1012,-1012,-1011,-1011,-1011,-1010,-1010,-1010,
                                 -1010,-1010,-1009,-1009,-1008,-1008,-1007,-1007,-1007,-1007,-1006,-1006,
                                 -1006,-1005,-1005,-1005,-1004,-1004,-1004,-1003,-1003,-1002,-1002,-1002,
                                 -1002,-1002,-1001,-1001,-1000, -999, -999, -999, -999, -999, -999, -999,
                                     -999, -998, -998, -997, -997, -997, -997, -996, -996, -996, -996, -996,
                                     -996, -995, -995, -994, -993, -993, -993, -993, -992, -992, -992, -992,
                                     -991, -991, -991, -990, -990, -989, -989, -988, -988, -988, -988, -987,
                                     -986, -986, -985, -985, -985, -985, -985, -985, -984, -984, -983, -983,
                                     -983, -982, -982, -981, -981, -981, -981, -980, -980, -980, -980, -980,
                                     -979, -979, -979, -979, -978, -978, -978, -977, -977, -976, -976, -976,
                                     -975, -975, -975, -974, -974, -974, -973, -973, -973, -973, -972, -972,
                                     -971, -971, -971, -970, -970, -970, -970, -969, -969, -969, -969, -968,
                                     -968, -968, -967, -967, -966, -966, -966, -965, -965, -965, -965, -965,
                                     -965, -965, -964, -964, -964, -964, -964, -963, -963, -963, -963, -962,
                                     -962, -962, -961, -961, -960, -960, -960, -960, -959, -959, -958, -958,
                                     -958, -957, -957, -957, -957, -956, -956, -955, -955, -955, -955, -955,
                                     -953, -953, -952, -952, -952, -951, -951, -951, -951, -951, -951, -950,
                                     -950, -950, -949, -949, -948, -948, -948, -948, -947, -947, -947, -946,
                                     -946, -946, -946, -945, -945, -945, -945, -944, -944, -944, -943, -943,
                                     -942, -942, -942, -942, -942, -941, -941, -940, -940, -940, -940, -940,
                                     -940, -939, -939, -939, -938, -937, -937, -937, -937, -936, -936, -936,
                                     -935, -935, -935, -934, -934, -934, -933, -933, -933, -933, -932, -932,
                                     -932, -931, -931, -930, -930, -930, -930, -929, -929, -929, -929, -929,
                                     -928, -928, -928, -927, -927, -927, -926, -926, -926, -926, -926, -926,
                                     -925, -925, -925, -925, -925, -925, -924, -924, -924, -923, -923, -923,
                                     -923, -923, -921, -921, -921, -921, -921, -920, -920, -920, -918, -918,
                                     -917, -917, -916, -916, -916, -916, -916, -916, -915, -915, -914, -914,
                                     -914, -914, -914, -914, -914, -913, -913, -912, -912, -912, -911, -911,
                                     -911, -910, -910, -910, -910, -910, -910, -910, -909, -909, -909, -909,
                                     -909, -909, -909, -908, -908, -908, -907, -907, -907, -907, -906, -906,
                                     -906, -906, -905, -905, -905, -905, -904, -904, -904, -903, -901, -901,
                                     -901, -900, -900, -899, -899, -899, -898, -898, -898, -898, -897, -897,
                                     -896, -896, -896, -896, -896, -896, -896, -896, -895, -895, -895, -895,
                                     -895, -894, -894, -894, -894, -894, -893, -893, -893, -893, -893, -892,
                                     -892, -891, -891, -891, -891, -891, -891, -891, -891, -890, -890, -890,
                                     -889, -889, -888, -888, -888, -888, -887, -887, -887, -886, -886, -886,
                                     -885, -885, -885, -885, -885, -884, -884, -883, -883, -882, -882, -882,
                                     -882, -881, -881, -881, -880, -880, -880, -880, -880, -879, -879, -879,
                                     -879, -878, -878, -878, -878, -878, -878, -878, -877, -877, -877, -877,
                                     -877, -877, -877, -876, -876, -876, -876, -876, -876, -876, -876, -875,
                                     -875, -875, -875, -874, -874, -874, -873, -873, -873, -873, -873, -872,
                                     -872, -871, -871, -871, -871, -871, -871, -871, -870, -870, -870, -870,
                                     -869, -869, -869, -869, -868, -868, -868, -868, -867, -867, -867, -867,
                                     -866, -866, -866, -866, -865, -865, -864, -864, -864, -864, -864, -864,
                                     -864, -864, -863, -863, -863, -862, -862, -862, -862, -861, -861, -861,
                                     -861, -861, -860, -860, -860, -859, -859, -859, -859, -859, -858, -858,
                                     -858, -858, -858, -857, -857, -857, -857, -856, -856, -855, -855, -855,
                                     -855, -855, -855, -855, -854, -854, -854, -854, -854, -853, -853, -853,
                                     -852, -851, -851, -851, -851, -851, -850, -850, -850, -850, -850, -850,
                                     -849, -849, -849, -848, -848, -848, -848, -847, -846, -846, -846, -846,
                                     -846, -846, -846, -846, -846, -846, -845, -845, -845, -844, -844, -844,
                                     -844, -844, -843, -843, -843, -843, -843, -842, -842, -842, -842, -842,
                                     -841, -841, -840, -840, -840, -840, -840, -839, -839, -839, -839, -839,
                                     -839, -837, -837, -837, -837, -837, -837, -837, -836, -836, -836, -836,
                                     -836, -836, -835, -835, -835, -835, -835, -835, -834, -834, -833, -833,
                                     -833, -833, -833, -832, -832, -832, -832, -831, -831, -831, -831, -831,
                                     -831, -830, -830, -830, -830, -830, -830, -830, -830, -829, -829, -829,
                                     -829, -828, -828, -828, -828, -828, -827, -827, -827, -827, -827, -827,
                                     -827, -827, -826, -826, -825, -825, -825, -825, -825, -824, -824, -824,
                                     -824, -824, -824, -823, -823, -823, -822, -822, -822, -822, -822, -822,
                                     -821, -821, -821, -820, -819, -819, -819, -819, -819, -819, -819, -818,
                                     -818, -818, -818, -817, -817, -817, -816, -816, -816, -816, -816, -816,
                                     -816, -816, -816, -815, -815, -815, -815, -815, -814, -814, -814, -813,
                                     -813, -813, -813, -813, -813, -813, -812, -812, -812, -811, -811, -811,
                                     -811, -811, -810, -810, -810, -810, -810, -809, -809, -809, -809, -809,
                                     -809, -809, -808, -808, -808, -807, -807, -807, -807, -807, -807, -807,
                                     -807, -807, -806, -806, -806, -805, -805, -805, -805, -804, -804, -804,
                                     -804, -803, -803, -802, -802, -802, -802, -802, -802, -802, -801, -801,
                                     -800, -800, -800, -800, -800, -799, -799, -799, -799, -799, -798, -798,
                                     -798, -798, -797, -797, -797, -797, -797, -796, -796, -796, -796, -796,
                                     -796, -796, -796, -796, -795, -795, -795, -795, -795, -795, -794, -794,
                                     -794, -794, -794, -794, -794, -793, -793, -792, -792, -792, -792, -792,
                                     -792, -792, -791, -791, -791, -791, -790, -790, -790, -790, -790, -790,
                                     -789, -789, -789, -789, -789, -788, -787, -787, -787, -787, -787, -787,
                                     -787, -786, -786, -786, -786, -785, -785, -785, -785, -784, -784, -784,
                                     -784, -784, -783, -783, -782, -782, -782, -782, -781, -781, -781, -781,
                                     -781, -781, -781, -781, -781, -780, -780, -780, -780, -780, -779, -779,
                                     -779, -779, -779, -779, -778, -778, -778, -778, -777, -777, -777, -777,
                                     -776, -776, -776, -776, -776, -776, -776, -776, -776, -775, -775, -775,
                                     -775, -774, -774, -774, -774, -773, -773, -773, -773, -773, -773, -773,
                                     -773, -773, -773, -773, -773, -773, -772, -771, -771, -771, -771, -770,
                                     -770, -770, -770, -770, -770, -770, -769, -769, -769, -769, -768, -768,
                                     -768, -768, -768, -768, -767, -767, -767, -767, -766, -766, -766, -766,
                                     -766, -766, -766, -765, -765, -765, -765, -764, -764, -764, -764, -764,
                                     -764, -763, -763, -763, -763, -763, -762, -762, -762, -762, -761, -761,
                                     -761, -761, -761, -761, -761, -760, -760, -760, -760, -760, -760, -760,
                                     -759, -759, -759, -758, -758, -758, -758, -758, -758, -757, -757, -757,
                                     -757, -757, -757, -756, -756, -756, -756, -756, -756, -756, -755, -755,
                                     -755, -755, -754, -754, -754, -754, -753, -753, -753, -753, -753, -753,
                                     -753, -752, -752, -752, -752, -751, -751, -751, -751, -751, -751, -751,
                                     -750, -750, -750, -749, -749, -749, -749, -748, -748, -748, -748, -748,
                                     -748, -748, -747, -747, -747, -746, -746, -745, -744, -744, -744, -744,
                                     -744, -744, -744, -743, -743, -743, -743, -743, -743, -743, -743, -742,
                                     -742, -741, -741, -741, -741, -741, -741, -741, -741, -741, -741, -740,
                                     -740, -740, -739, -739, -739, -739, -739, -739, -739, -739, -739, -739,
                                     -738, -738, -738, -738, -738, -738, -738, -737, -737, -736, -736, -736,
                                     -736, -736, -736, -736, -736, -736, -736, -735, -735, -735, -735, -735,
                                     -734, -734, -734, -734, -734, -733, -733, -733, -733, -733, -732, -732,
                                     -732, -732, -732, -731, -731, -731, -730, -730, -730, -730, -730, -730,
                                     -729, -729, -729, -728, -728, -727, -727, -727, -727, -727, -726, -726,
                                     -726, -726, -726, -726, -726, -726, -726, -725, -725, -725, -725, -724,
                                     -724, -724, -724, -724, -724, -723, -723, -723, -722, -722, -722, -722,
                                     -722, -722, -722, -722, -722, -721, -721, -721, -721, -721, -721, -721,
                                     -721, -721, -720, -720, -720, -720, -719, -719, -719, -719, -719, -719,
                                     -719, -719, -719, -719, -718, -718, -718, -718, -717, -717, -717, -717,
                                     -716, -716, -716, -716, -715, -715, -715, -715, -714, -714, -714, -714,
                                     -714, -713, -713, -713, -713, -713, -712, -712, -712, -712, -712, -712,
                                     -712, -711, -711, -711, -711, -711, -710, -710, -710, -710, -710, -710,
                                     -710, -710, -709, -709, -709, -709, -709, -709, -708, -708, -708, -708,
                                     -708, -707, -707, -707, -707, -707, -707, -707, -707, -706, -706, -706,
                                     -706, -706, -706, -705, -705, -705, -705, -705, -705, -705, -705, -704,
                                     -704, -704, -704, -704, -704, -704, -704, -703, -703, -703, -703, -702,
                                     -702, -702, -702, -702, -702, -702, -702, -702, -702, -702, -701, -701,
                                     -701, -701, -701, -700, -700, -699, -699, -699, -699, -699, -699, -699,
                                     -699, -698, -698, -698, -697, -697, -697, -697, -697, -697, -697, -697,
                                     -696, -696, -696, -696, -696, -696, -695, -695, -695, -694, -694, -694,
                                     -694, -694, -694, -694, -694, -694, -694, -694, -694, -694, -694, -693,
                                     -693, -693, -693, -692, -692, -692, -692, -692, -692, -692, -692, -691,
                                     -691, -691, -691, -691, -690, -690, -690, -690, -690, -689, -689, -689,
                                     -689, -689, -689, -689, -689, -688, -688, -688, -688, -688, -688, -687,
                                     -687, -687, -687, -687, -687, -687, -686, -686, -686, -685, -685, -685,
                                     -685, -685, -685, -684, -684, -684, -684, -684, -684, -683, -683, -683,
                                     -683, -683, -682, -682, -682, -682, -682, -682, -682, -681, -681, -681,
                                     -681, -681, -680, -680, -680, -680, -679, -679, -679, -679, -679, -679,
                                     -678, -678, -678, -678, -678, -678, -678, -678, -678, -677, -677, -677,
                                     -677, -677, -677, -676, -676, -676, -676, -676, -676, -675, -675, -675,
                                     -675, -675, -675, -675, -675, -674, -674, -674, -674, -674, -674, -674,
                                     -674, -673, -673, -673, -673, -673, -673, -673, -673, -672, -672, -672,
                                     -672, -672, -672, -672, -671, -671, -670, -670, -670, -670, -670, -670,
                                     -670, -669, -669, -669, -669, -669, -668, -668, -667, -667, -667, -667,
                                     -667, -667, -667, -667, -667, -666, -666, -666, -666, -666, -666, -665,
                                     -665, -665, -665, -665, -665, -664, -664, -664, -664, -664, -664, -663,
                                     -663, -663, -663, -663, -663, -663, -663, -662, -662, -662, -662, -662,
                                     -662, -662, -661, -661, -661, -661, -661, -661, -661, -661, -660, -660,
                                     -660, -660, -660, -660, -660, -660, -660, -659, -659, -659, -659, -659,
                                     -658, -658, -658, -658, -658, -658, -658, -657, -657, -657, -657, -657,
                                     -657, -657, -657, -656, -656, -656, -656, -656, -655, -655, -655, -655,
                                     -655, -655, -655, -655, -654, -654, -654, -654, -654, -654, -654, -654,
                                     -653, -653, -653, -653, -653, -653, -652, -652, -652, -652, -652, -651,
                                     -651, -651, -651, -650, -650, -650, -650, -650, -650, -650, -649, -649,
                                     -649, -649, -649, -649, -649, -648, -648, -648, -648, -648, -648, -648,
                                     -647, -647, -647, -647, -647, -647, -646, -646, -646, -646, -646, -646,
                                     -646, -645, -645, -645, -645, -645, -645, -645, -644, -644, -644, -644,
                                     -644, -644, -644, -644, -644, -643, -643, -643, -643, -643, -643, -643,
                                     -642, -642, -642, -642, -641, -641, -641, -641, -641, -641, -641, -641,
                                     -640, -640, -640, -640, -640, -640, -639, -639, -639, -639, -639, -638,
                                     -638, -638, -638, -638, -638, -638, -638, -638, -638, -638, -638, -637,
                                     -637, -637, -637, -637, -637, -637, -636, -636, -636, -636, -636, -636,
                                     -636, -635, -635, -635, -635, -634, -634, -634, -634, -633, -633, -633,
                                     -633, -633, -633, -633, -633, -632, -632, -632, -631, -631, -631, -631,
                                     -631, -631, -631, -630, -630, -630, -630, -630, -630, -630, -630, -630,
                                     -630, -630, -629, -629, -629, -629, -629, -629, -629, -628, -628, -628,
                                     -628, -627, -627, -627, -627, -627, -627, -626, -626, -626, -625, -625,
                                     -625, -625, -625, -625, -624, -624, -624, -624, -624, -624, -624, -624,
                                     -624, -624, -623, -623, -623, -623, -623, -623, -623, -623, -623, -622,
                                     -622, -622, -622, -622, -622, -621, -621, -621, -621, -621, -620, -620,
                                     -620, -620, -620, -619, -619, -619, -619, -619, -619, -619, -619, -618,
                                     -618, -618, -618, -618, -618, -618, -618, -617, -617, -617, -617, -616,
                                     -616, -616, -616, -615, -615, -615, -615, -615, -614, -614, -614, -614,
                                     -614, -614, -613, -613, -613, -613, -613, -613, -613, -613, -612, -612,
                                     -612, -612, -612, -612, -612, -612, -612, -612, -612, -612, -611, -611,
                                     -611, -611, -611, -611, -610, -610, -610, -610, -610, -610, -609, -609,
                                     -609, -609, -609, -609, -609, -608, -608, -608, -608, -608, -608, -608,
                                     -607, -607, -607, -607, -607, -607, -607, -606, -606, -606, -606, -605,
                                     -605, -605, -605, -605, -605, -604, -604, -604, -604, -604, -604, -604,
                                     -603, -603, -603, -603, -603, -603, -602, -602, -602, -602, -602, -602,
                                     -602, -601, -601, -601, -601, -601, -601, -601, -601, -601, -600, -600,
                                     -600, -600, -600, -599, -599, -599, -599, -599, -599, -599, -599, -598,
                                     -598, -598, -598, -598, -597, -597, -597, -597, -597, -597, -597, -597,
                                     -597, -597, -596, -596, -596, -596, -596, -596, -596, -596, -596, -596,
                                     -596, -595, -595, -595, -595, -594, -594, -594, -594, -594, -593, -593,
                                     -593, -593, -593, -593, -593, -593, -593, -593, -593, -592, -592, -592,
                                     -592, -591, -591, -591, -591, -591, -591, -591, -591, -591, -590, -590,
                                     -590, -590, -590, -590, -590, -589, -589, -589, -589, -589, -589, -589,
                                     -588, -588, -588, -588, -588, -588, -587, -587, -587, -587, -586, -586,
                                     -586, -586, -586, -586, -586, -586, -585, -585, -585, -585, -585, -584,
                                     -584, -584, -584, -584, -583, -583, -583, -583, -582, -582, -582, -582,
                                     -582, -582, -582, -582, -582, -581, -581, -581, -581, -581, -581, -581,
                                     -580, -580, -580, -580, -580, -580, -580, -580, -579, -579, -579, -579,
                                     -579, -579, -579, -579, -579, -579, -579, -579, -579, -578, -578, -578,
                                     -578, -578, -578, -577, -577, -577, -577, -577, -576, -576, -576, -576,
                                     -576, -576, -576, -575, -575, -575, -575, -575, -575, -574, -574, -574,
                                     -574, -574, -573, -573, -573, -573, -573, -573, -573, -573, -573, -572,
                                     -572, -572, -572, -572, -572, -572, -572, -572, -571, -571, -571, -571,
                                     -571, -571, -571, -570, -570, -570, -570, -570, -570, -569, -569, -569,
                                     -569, -569, -569, -569, -568, -568, -568, -568, -568, -567, -567, -567,
                                     -567, -567, -567, -567, -566, -566, -566, -566, -566, -566, -565, -565,
                                     -565, -565, -565, -564, -564, -564, -564, -564, -564, -564, -564, -563,
                                     -563, -563, -563, -563, -563, -563, -562, -562, -562, -562, -562, -562,
                                     -562, -562, -562, -561, -561, -561, -561, -561, -561, -561, -561, -561,
                                     -560, -560, -560, -560, -560, -560, -559, -559, -559, -559, -559, -559,
                                     -558, -558, -558, -558, -558, -558, -558, -558, -557, -557, -557, -557,
                                     -557, -557, -557, -557, -556, -556, -556, -556, -556, -556, -556, -555,
                                     -555, -555, -555, -555, -555, -555, -555, -555, -554, -554, -554, -554,
                                     -554, -554, -554, -553, -553, -553, -553, -553, -553, -553, -553, -552,
                                     -552, -552, -552, -552, -552, -552, -552, -552, -552, -551, -551, -550,
                                     -550, -550, -550, -550, -549, -549, -549, -549, -549, -549, -549, -549,
                                     -549, -549, -548, -548, -548, -548, -548, -548, -548, -548, -547, -547,
                                     -547, -547, -546, -546, -546, -546, -546, -546, -546, -545, -545, -545,
                                     -545, -545, -545, -544, -544, -544, -544, -544, -544, -544, -544, -544,
                                     -543, -543, -543, -543, -543, -543, -542, -542, -542, -542, -542, -541,
                                     -541, -541, -541, -541, -541, -541, -541, -540, -540, -540, -540, -540,
                                     -540, -540, -540, -540, -540, -539, -539, -539, -539, -539, -539, -539,
                                     -539, -539, -539, -539, -538, -538, -538, -538, -538, -538, -538, -538,
                                     -538, -538, -538, -538, -537, -537, -537, -537, -537, -537, -537, -537,
                                     -536, -536, -536, -536, -536, -536, -536, -536, -536, -536, -535, -535,
                                     -535, -535, -535, -535, -535, -535, -535, -534, -534, -534, -534, -533,
                                     -533, -533, -533, -532, -532, -532, -532, -532, -532, -531, -531, -531,
                                     -531, -531, -531, -531, -531, -530, -530, -530, -530, -530, -530, -530,
                                     -530, -530, -530, -529, -529, -529, -528, -528, -528, -528, -528, -528,
                                     -528, -528, -528, -528, -528, -527, -527, -527, -527, -527, -527, -527,
                                     -527, -527, -527, -526, -526, -526, -526, -526, -525, -525, -525, -525,
                                     -525, -525, -525, -525, -524, -524, -524, -524, -524, -524, -524, -523,
                                     -523, -523, -523, -523, -523, -523, -522, -522, -522, -522, -522, -522,
                                     -522, -522, -522, -522, -521, -521, -521, -521, -521, -521, -521, -521,
                                     -521, -520, -520, -520, -520, -520, -520, -520, -520, -520, -520, -520,
                                     -520, -520, -519, -519, -519, -519, -519, -518, -518, -518, -518, -518,
                                     -518, -518, -518, -517, -517, -517, -517, -517, -517, -517, -517, -517,
                                     -517, -516, -516, -516, -515, -515, -515, -515, -515, -515, -515, -515,
                                     -515, -514, -514, -514, -514, -513, -513, -513, -513, -513, -513, -513,
                                     -512, -512, -512, -512, -512, -512, -511, -511, -511, -511, -511, -511,
                                     -511, -511, -510, -510, -510, -510, -510, -510, -510, -509, -509, -509,
                                     -509, -509, -509, -509, -509, -509, -509, -508, -508, -508, -508, -508,
                                     -508, -508, -508, -508, -507, -507, -507, -507, -507, -506, -506, -506,
                                     -506, -506, -506, -506, -506, -506, -505, -505, -505, -505, -505, -505,
                                     -505, -505, -505, -504, -504, -504, -504, -504, -504, -504, -504, -504,
                                     -503, -503, -503, -503, -503, -503, -503, -503, -503, -502, -502, -502,
                                     -502, -502, -502, -502, -502, -502, -502, -502, -501, -501, -501, -501,
                                     -501, -501, -501, -501, -500, -500, -500, -500, -500, -500, -500, -500,
                                     -500, -499, -499, -499, -499, -499, -499, -498, -498, -498, -498, -498,
                                     -498, -498, -497, -497, -497, -497, -497, -497, -497, -497, -497, -497,
                                     -496, -496, -496, -496, -496, -496, -496, -496, -496, -495, -495, -495,
                                     -495, -495, -495, -495, -495, -495, -494, -494, -494, -493, -493, -493,
                                     -493, -493, -493, -493, -493, -492, -492, -492, -492, -492, -491, -491,
                                     -491, -491, -491, -491, -491, -491, -491, -491, -491, -491, -490, -490,
                                     -490, -490, -490, -490, -490, -490, -490, -490, -490, -489, -489, -489,
                                     -489, -489, -489, -489, -488, -488, -488, -488, -488, -488, -488, -488,
                                     -488, -488, -488, -488, -488, -487, -487, -487, -487, -486, -486, -486,
                                     -486, -486, -486, -486, -485, -485, -485, -485, -485, -485, -485, -485,
                                     -485, -484, -484, -484, -484, -484, -484, -483, -483, -483, -483, -483,
                                     -483, -482, -482, -482, -482, -482, -482, -482, -482, -482, -481, -481,
                                     -481, -481, -481, -481, -481, -481, -481, -481, -480, -480, -480, -480,
                                     -480, -480, -480, -480, -479, -479, -479, -479, -479, -479, -479, -479,
                                     -479, -479, -479, -478, -478, -478, -478, -478, -478, -478, -478, -478,
                                     -478, -478, -477, -477, -477, -477, -476, -476, -476, -476, -476, -476,
                                     -476, -476, -476, -476, -476, -475, -475, -475, -475, -475, -475, -475,
                                     -475, -475, -475, -474, -474, -474, -474, -473, -473, -473, -473, -473,
                                     -473, -473, -472, -472, -472, -472, -472, -472, -472, -471, -471, -471,
                                     -471, -471, -471, -471, -471, -471, -471, -471, -471, -470, -470, -470,
                                     -470, -470, -469, -469, -469, -469, -469, -469, -469, -468, -468, -468,
                                     -468, -468, -468, -468, -468, -468, -468, -468, -467, -467, -467, -466,
                                     -466, -466, -466, -466, -466, -466, -466, -465, -465, -465, -465, -465,
                                     -464, -464, -464, -464, -464, -464, -464, -464, -464, -464, -464, -464,
                                     -463, -463, -463, -463, -463, -463, -463, -463, -463, -463, -462, -462,
                                     -462, -462, -462, -462, -462, -462, -462, -461, -461, -461, -461, -461,
                                     -461, -461, -461, -461, -461, -461, -461, -461, -461, -461, -461, -460,
                                     -460, -460, -460, -460, -460, -460, -460, -460, -459, -459, -459, -459,
                                     -459, -458, -458, -458, -458, -458, -458, -458, -458, -458, -457, -457,
                                     -457, -457, -457, -457, -456, -456, -456, -456, -455, -455, -455, -455,
                                     -455, -454, -454, -454, -454, -454, -454, -454, -453, -453, -453, -453,
                                     -453, -453, -453, -453, -453, -453, -452, -452, -452, -452, -452, -452,
                                     -452, -452, -452, -451, -451, -451, -451, -451, -451, -450, -450, -450,
                                     -450, -450, -450, -450, -450, -450, -449, -449, -449, -449, -449, -449,
                                     -449, -449, -449, -449, -448, -448, -448, -448, -448, -448, -447, -447,
                                     -447, -447, -447, -447, -447, -446, -446, -446, -446, -446, -446, -446,
                                     -446, -446, -446, -446, -446, -446, -445, -445, -445, -445, -445, -445,
                                     -445, -445, -445, -445, -445, -444, -444, -444, -444, -444, -444, -444,
                                     -444, -443, -443, -443, -443, -443, -443, -443, -443, -443, -443, -443,
                                     -442, -442, -442, -442, -442, -442, -442, -441, -441, -441, -441, -441,
                                     -441, -441, -441, -441, -441, -441, -440, -440, -440, -440, -440, -440,
                                     -440, -439, -439, -439, -438, -438, -438, -438, -438, -438, -438, -438,
                                     -438, -438, -437, -437, -437, -437, -437, -437, -437, -437, -436, -436,
                                     -436, -436, -436, -436, -436, -436, -436, -436, -435, -435, -435, -435,
                                     -435, -435, -435, -435, -434, -434, -434, -434, -434, -434, -433, -433,
                                     -433, -433, -433, -433, -433, -433, -433, -433, -433, -432, -432, -432,
                                     -432, -432, -432, -432, -432, -432, -431, -431, -431, -431, -431, -431,
                                     -431, -431, -431, -431, -431, -431, -430, -430, -430, -430, -430, -430,
                                     -430, -430, -430, -429, -429, -429, -429, -429, -429, -428, -428, -428,
                                     -428, -428, -428, -428, -428, -428, -428, -427, -427, -427, -427, -427,
                                     -427, -427, -427, -427, -427, -427, -426, -426, -426, -426, -426, -426,
                                     -425, -425, -425, -425, -425, -425, -425, -425, -425, -425, -424, -424,
                                     -424, -424, -424, -424, -424, -424, -423, -423, -423, -423, -423, -423,
                                     -423, -423, -422, -422, -422, -422, -422, -421, -421, -421, -421, -421,
                                     -421, -421, -421, -421, -421, -421, -420, -420, -420, -420, -420, -420,
                                     -420, -420, -420, -420, -420, -420, -419, -419, -419, -419, -419, -419,
                                     -419, -418, -418, -418, -418, -417, -417, -417, -417, -417, -417, -417,
                                     -416, -416, -416, -416, -416, -416, -416, -416, -416, -415, -415, -415,
                                     -415, -415, -415, -415, -415, -415, -415, -414, -414, -414, -414, -414,
                                     -414, -414, -413, -413, -413, -413, -413, -413, -413, -413, -413, -412,
                                     -412, -412, -412, -412, -412, -412, -412, -412, -412, -412, -412, -411,
                                     -411, -411, -411, -411, -411, -411, -410, -410, -410, -410, -410, -410,
                                     -410, -410, -410, -409, -409, -409, -409, -409, -408, -408, -408, -408,
                                     -408, -408, -408, -408, -407, -407, -407, -407, -407, -407, -406, -406,
                                     -406, -406, -406, -406, -406, -406, -406, -406, -406, -406, -405, -405,
                                     -405, -405, -405, -405, -405, -405, -405, -405, -405, -404, -404, -404,
                                     -404, -404, -404, -404, -404, -403, -403, -403, -403, -403, -403, -403,
                                     -403, -403, -403, -403, -403, -403, -403, -403, -403, -403, -403, -403,
                                     -403, -403, -402, -402, -402, -402, -402, -402, -402, -401, -401, -401,
                                     -401, -401, -401, -401, -401, -400, -400, -400, -400, -400, -400, -399,
                                     -399, -399, -398, -398, -398, -398, -398, -398, -398, -398, -398, -398,
                                     -397, -397, -397, -397, -397, -397, -397, -397, -397, -397, -396, -396,
                                     -396, -396, -396, -396, -396, -396, -395, -395, -395, -395, -395, -394,
                                     -394, -394, -394, -394, -394, -394, -394, -394, -394, -394, -394, -394,
                                     -394, -393, -393, -393, -393, -393, -393, -393, -393, -392, -392, -392,
                                     -392, -392, -392, -392, -392, -392, -391, -391, -391, -391, -391, -391,
                                     -391, -391, -391, -390, -390, -390, -390, -390, -390, -389, -389, -389,
                                     -389, -389, -389, -389, -389, -389, -389, -389, -389, -388, -388, -388,
                                     -388, -388, -388, -388, -388, -388, -388, -387, -387, -387, -387, -387,
                                     -387, -387, -387, -387, -387, -386, -386, -386, -386, -386, -386, -386,
                                     -386, -386, -386, -386, -386, -385, -385, -385, -385, -385, -385, -385,
                                     -385, -385, -385, -385, -384, -384, -384, -384, -384, -384, -384, -384,
                                     -384, -384, -383, -383, -383, -383, -383, -383, -383, -383, -383, -383,
                                     -382, -382, -382, -382, -382, -382, -382, -382, -382, -382, -382, -381,
                                     -381, -381, -381, -381, -380, -380, -380, -380, -380, -380, -379, -379,
                                     -379, -379, -379, -378, -378, -378, -378, -378, -378, -378, -378, -378,
                                     -378, -378, -378, -378, -377, -377, -377, -377, -377, -376, -376, -376,
                                     -376, -376, -376, -376, -376, -376, -376, -376, -376, -375, -375, -375,
                                     -375, -375, -375, -375, -374, -374, -374, -374, -374, -374, -374, -374,
                                     -373, -373, -373, -373, -373, -373, -373, -373, -373, -372, -372, -372,
                                     -372, -372, -372, -372, -372, -372, -372, -372, -372, -371, -371, -371,
                                     -371, -371, -371, -371, -371, -371, -371, -370, -370, -370, -370, -370,
                                     -370, -370, -370, -369, -369, -369, -369, -369, -369, -369, -369, -369,
                                     -369, -369, -369, -369, -369, -368, -368, -368, -368, -368, -368, -368,
                                     -368, -368, -368, -368, -368, -368, -367, -367, -367, -367, -367, -367,
                                     -367, -367, -367, -367, -367, -367, -367, -366, -366, -366, -366, -366,
                                     -366, -366, -366, -366, -366, -366, -366, -365, -365, -365, -365, -365,
                                     -365, -364, -364, -364, -364, -364, -364, -364, -363, -363, -363, -363,
                                     -363, -363, -363, -363, -363, -363, -362, -362, -362, -362, -362, -361,
                                     -361, -361, -361, -361, -361, -361, -361, -361, -361, -361, -361, -361,
                                     -361, -360, -360, -360, -360, -360, -359, -359, -359, -359, -359, -359,
                                     -358, -358, -358, -358, -358, -358, -358, -358, -358, -358, -358, -357,
                                     -357, -357, -357, -357, -357, -356, -356, -356, -356, -356, -356, -356,
                                     -356, -356, -356, -356, -355, -355, -355, -355, -355, -355, -355, -355,
                                     -355, -355, -355, -355, -355, -355, -355, -354, -354, -354, -354, -354,
                                     -354, -354, -354, -354, -354, -353, -353, -353, -353, -353, -353, -353,
                                     -353, -353, -352, -352, -352, -352, -352, -352, -352, -352, -351, -351,
                                     -351, -351, -351, -351, -351, -351, -351, -351, -351, -350, -350, -350,
                                     -350, -350, -350, -350, -350, -350, -349, -349, -349, -349, -349, -349,
                                     -349, -349, -349, -349, -349, -349, -348, -348, -348, -348, -348, -348,
                                     -348, -347, -347, -347, -347, -347, -347, -347, -347, -347, -347, -346,
                                     -346, -346, -346, -346, -346, -346, -346, -346, -345, -345, -345, -345,
                                     -345, -345, -345, -345, -344, -344, -344, -344, -344, -344, -344, -344,
                                     -343, -343, -343, -343, -343, -343, -343, -342, -342, -342, -342, -342,
                                     -341, -341, -341, -341, -341, -341, -341, -341, -341, -341, -341, -340,
                                     -340, -340, -340, -340, -340, -340, -340, -340, -340, -340, -340, -339,
                                     -339, -339, -339, -339, -339, -339, -339, -338, -338, -338, -338, -338,
                                     -338, -338, -338, -338, -338, -337, -337, -337, -337, -337, -337, -337,
                                     -337, -337, -337, -336, -336, -336, -336, -336, -336, -336, -336, -336,
                                     -336, -336, -335, -335, -335, -335, -335, -335, -335, -335, -335, -335,
                                     -334, -334, -334, -334, -334, -334, -334, -334, -334, -334, -334, -334,
                                     -333, -333, -333, -333, -333, -333, -333, -333, -333, -333, -333, -333,
                                     -332, -332, -332, -332, -332, -332, -332, -332, -332, -332, -332, -332,
                                     -332, -331, -331, -331, -331, -331, -331, -331, -331, -330, -330, -330,
                                     -330, -330, -330, -330, -329, -329, -329, -329, -329, -329, -329, -328,
                                     -328, -328, -328, -328, -328, -328, -328, -328, -328, -328, -328, -328,
                                     -328, -327, -327, -327, -327, -327, -327, -327, -327, -326, -326, -326,
                                     -326, -326, -326, -326, -326, -326, -326, -326, -325, -325, -325, -325,
                                     -325, -325, -325, -325, -324, -324, -324, -324, -324, -324, -324, -324,
                                     -324, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -322,
                                     -322, -322, -322, -322, -322, -322, -322, -322, -321, -321, -321, -321,
                                     -321, -321, -321, -321, -321, -321, -321, -321, -320, -320, -320, -320,
                                     -320, -320, -319, -319, -319, -319, -319, -319, -319, -319, -319, -319,
                                     -319, -318, -318, -318, -318, -318, -318, -318, -318, -318, -318, -317,
                                     -317, -317, -317, -317, -317, -317, -317, -316, -316, -316, -316, -316,
                                     -316, -316, -316, -316, -316, -316, -316, -316, -316, -316, -315, -315,
                                     -315, -315, -315, -315, -315, -315, -315, -315, -315, -315, -314, -314,
                                     -314, -314, -314, -314, -314, -314, -314, -313, -313, -313, -313, -313,
                                     -313, -313, -312, -312, -312, -312, -312, -312, -312, -312, -312, -311,
                                     -311, -311, -311, -311, -310, -310, -310, -310, -310, -310, -310, -310,
                                     -310, -310, -309, -309, -309, -309, -309, -309, -309, -309, -309, -308,
                                     -308, -308, -308, -308, -308, -308, -308, -308, -308, -308, -307, -307,
                                     -307, -307, -307, -307, -307, -307, -307, -307, -307, -307, -307, -306,
                                     -306, -306, -306, -306, -306, -306, -306, -305, -305, -305, -305, -305,
                                     -305, -305, -305, -305, -305, -305, -305, -305, -305, -305, -304, -304,
                                     -304, -304, -304, -303, -303, -303, -303, -303, -303, -303, -303, -303,
                                     -303, -303, -302, -302, -302, -302, -302, -302, -302, -302, -302, -302,
                                     -301, -301, -301, -301, -301, -301, -301, -301, -300, -300, -300, -300,
                                     -300, -300, -300, -300, -299, -299, -299, -299, -299, -299, -298, -298,
                                     -298, -298, -298, -298, -298, -298, -298, -297, -297, -297, -297, -297,
                                     -297, -297, -297, -297, -297, -296, -296, -296, -296, -296, -296, -296,
                                     -296, -296, -296, -296, -295, -295, -295, -295, -295, -295, -295, -295,
                                     -295, -295, -295, -295, -294, -294, -294, -294, -294, -294, -294, -294,
                                     -294, -294, -294, -294, -294, -293, -293, -293, -293, -293, -293, -293,
                                     -293, -292, -292, -292, -292, -292, -292, -292, -292, -291, -291, -291,
                                     -291, -291, -291, -291, -291, -291, -291, -291, -291, -291, -291, -291,
                                     -291, -291, -291, -291, -290, -290, -290, -290, -290, -289, -289, -289,
                                     -289, -289, -289, -289, -289, -289, -288, -288, -288, -288, -288, -288,
                                     -288, -288, -288, -288, -288, -288, -288, -287, -287, -287, -287, -287,
                                     -287, -286, -286, -286, -286, -286, -286, -286, -286, -286, -286, -285,
                                     -285, -285, -285, -285, -285, -285, -285, -285, -285, -285, -285, -285,
                                     -285, -285, -284, -284, -284, -284, -284, -284, -283, -283, -283, -283,
                                     -283, -283, -283, -283, -283, -282, -282, -282, -282, -282, -282, -281,
                                     -281, -281, -281, -281, -281, -281, -280, -280, -280, -280, -280, -280,
                                     -280, -280, -280, -279, -279, -279, -279, -279, -279, -278, -278, -278,
                                     -278, -278, -278, -278, -278, -278, -278, -278, -278, -278, -278, -277,
                                     -277, -277, -277, -277, -277, -277, -277, -277, -276, -276, -276, -276,
                                     -276, -276, -276, -276, -276, -276, -276, -275, -275, -275, -275, -275,
                                     -275, -275, -275, -275, -275, -275, -275, -275, -275, -274, -274, -274,
                                     -274, -274, -274, -274, -274, -273, -273, -273, -273, -273, -273, -273,
                                     -273, -273, -273, -273, -273, -273, -272, -272, -272, -272, -272, -272,
                                     -272, -272, -272, -272, -272, -272, -272, -271, -271, -271, -271, -271,
                                     -271, -271, -271, -271, -271, -271, -270, -270, -270, -270, -270, -270,
                                     -270, -270, -270, -269, -269, -269, -269, -269, -269, -269, -269, -268,
                                     -268, -268, -268, -268, -268, -268, -268, -268, -268, -267, -267, -267,
                                     -267, -267, -267, -267, -267, -266, -266, -266, -266, -266, -266, -266,
                                     -266, -266, -266, -266, -266, -266, -266, -265, -265, -265, -265, -265,
                                     -265, -264, -264, -264, -264, -264, -264, -264, -264, -264, -263, -263,
                                     -263, -263, -263, -263, -263, -263, -263, -263, -263, -262, -262, -262,
                                     -262, -262, -262, -262, -262, -261, -261, -261, -261, -261, -261, -261,
                                     -261, -261, -261, -261, -260, -260, -260, -260, -260, -260, -260, -260,
                                     -260, -260, -259, -259, -259, -259, -259, -258, -258, -258, -258, -258,
                                     -258, -258, -258, -258, -258, -258, -257, -257, -257, -257, -257, -257,
                                     -257, -257, -257, -257, -257, -256, -256, -256, -256, -256, -256, -256,
                                     -256, -256, -256, -256, -256, -256, -255, -255, -255, -255, -255, -255,
                                     -255, -255, -255, -255, -255, -255, -255, -254, -254, -254, -254, -254,
                                     -254, -254, -254, -253, -253, -253, -253, -253, -253, -253, -253, -253,
                                     -253, -253, -252, -252, -252, -252, -252, -252, -252, -252, -251, -251,
                                     -251, -251, -251, -251, -251, -251, -251, -251, -251, -251, -251, -251,
                                     -251, -250, -250, -250, -250, -250, -249, -249, -249, -249, -249, -249,
                                     -248, -248, -248, -248, -248, -248, -248, -248, -247, -247, -247, -247,
                                     -247, -247, -247, -247, -247, -247, -247, -247, -247, -246, -246, -246,
                                     -246, -246, -246, -246, -246, -246, -246, -246, -245, -245, -245, -245,
                                     -245, -245, -245, -245, -245, -244, -244, -244, -244, -244, -244, -244,
                                     -244, -244, -244, -243, -243, -243, -243, -243, -243, -243, -243, -242,
                                     -242, -242, -242, -242, -242, -242, -242, -241, -241, -241, -241, -241,
                                     -241, -241, -241, -241, -241, -241, -240, -240, -240, -240, -240, -240,
                                     -240, -239, -239, -239, -239, -239, -239, -239, -239, -239, -238, -238,
                                     -238, -238, -238, -238, -238, -237, -237, -237, -237, -237, -237, -237,
                                     -237, -236, -236, -236, -236, -236, -236, -236, -236, -236, -236, -236,
                                     -236, -235, -235, -235, -235, -235, -235, -235, -235, -235, -235, -234,
                                     -234, -234, -234, -234, -234, -234, -234, -234, -234, -234, -234, -234,
                                     -234, -234, -234, -233, -233, -233, -233, -233, -233, -233, -233, -233,
                                     -233, -233, -232, -232, -232, -232, -232, -232, -232, -232, -232, -231,
                                     -231, -231, -231, -231, -231, -231, -231, -231, -231, -231, -230, -230,
                                     -230, -230, -230, -230, -230, -230, -230, -230, -230, -230, -230, -230,
                                     -230, -229, -229, -229, -229, -229, -229, -229, -229, -229, -228, -228,
                                     -228, -228, -228, -228, -228, -228, -228, -228, -228, -227, -227, -227,
                                     -227, -227, -227, -226, -226, -226, -226, -226, -226, -226, -226, -226,
                                     -226, -226, -226, -226, -226, -225, -225, -225, -225, -225, -225, -225,
                                     -225, -224, -224, -224, -224, -224, -224, -224, -224, -224, -223, -223,
                                     -223, -223, -223, -223, -223, -223, -223, -223, -223, -222, -222, -222,
                                     -222, -222, -221, -221, -221, -221, -221, -221, -221, -221, -221, -221,
                                     -221, -221, -221, -220, -220, -220, -220, -220, -220, -220, -220, -220,
                                     -219, -219, -219, -219, -219, -219, -219, -218, -218, -218, -218, -218,
                                     -218, -218, -218, -218, -218, -218, -218, -218, -217, -217, -217, -217,
                                     -217, -217, -217, -217, -217, -216, -216, -216, -216, -216, -216, -216,
                                     -216, -216, -216, -216, -216, -216, -215, -215, -215, -215, -215, -215,
                                     -215, -215, -215, -215, -215, -215, -214, -214, -214, -214, -214, -214,
                                     -214, -214, -214, -213, -213, -213, -213, -213, -213, -213, -213, -213,
                                     -213, -213, -212, -212, -212, -212, -212, -212, -212, -212, -212, -212,
                                     -212, -212, -211, -211, -211, -211, -211, -211, -211, -211, -210, -210,
                                     -210, -210, -210, -210, -210, -210, -210, -209, -209, -209, -209, -209,
                                     -209, -209, -209, -208, -208, -208, -208, -208, -208, -208, -208, -208,
                                     -208, -207, -207, -207, -207, -207, -207, -207, -207, -207, -207, -207,
                                     -207, -206, -206, -206, -206, -206, -206, -206, -206, -206, -206, -206,
                                     -206, -205, -205, -205, -205, -205, -205, -205, -205, -205, -204, -204,
                                     -204, -204, -204, -204, -204, -204, -204, -204, -204, -203, -203, -203,
                                     -203, -203, -203, -203, -203, -203, -203, -203, -203, -202, -202, -202,
                                     -202, -202, -202, -202, -202, -202, -202, -201, -201, -201, -201, -201,
                                     -201, -201, -201, -201, -201, -201, -200, -200, -200, -200, -200, -200,
                                     -200, -200, -200, -200, -199, -199, -199, -199, -199, -199, -199, -199,
                                     -199, -199, -199, -199, -199, -199, -199, -199, -198, -198, -198, -198,
                                     -198, -198, -198, -197, -197, -197, -197, -197, -197, -197, -197, -197,
                                     -197, -197, -197, -197, -197, -196, -196, -196, -196, -196, -196, -195,
                                     -195, -195, -195, -195, -195, -195, -195, -195, -195, -195, -195, -194,
                                     -194, -194, -194, -194, -194, -194, -194, -194, -194, -194, -193, -193,
                                     -193, -193, -193, -193, -193, -193, -192, -192, -192, -192, -192, -192,
                                     -192, -192, -192, -191, -191, -191, -191, -191, -191, -191, -191, -191,
                                     -190, -190, -190, -190, -190, -190, -190, -190, -190, -190, -189, -189,
                                     -189, -189, -189, -189, -189, -189, -189, -189, -188, -188, -188, -188,
                                     -188, -188, -188, -188, -188, -188, -188, -187, -187, -187, -187, -187,
                                     -187, -187, -187, -187, -187, -187, -187, -186, -186, -186, -186, -186,
                                     -186, -186, -186, -186, -186, -186, -186, -186, -185, -185, -185, -185,
                                     -185, -185, -185, -185, -185, -185, -185, -185, -184, -184, -184, -184,
                                     -184, -184, -184, -184, -184, -184, -184, -184, -184, -184, -184, -184,
                                     -183, -183, -183, -183, -183, -183, -183, -182, -182, -182, -182, -182,
                                     -182, -182, -182, -182, -182, -181, -181, -181, -181, -181, -181, -181,
                                     -181, -180, -180, -180, -180, -180, -180, -180, -180, -180, -180, -180,
                                     -179, -179, -179, -179, -179, -179, -179, -179, -179, -178, -178, -178,
                                     -178, -178, -178, -178, -177, -177, -177, -177, -177, -177, -177, -177,
                                     -177, -177, -177, -177, -177, -177, -177, -177, -176, -176, -176, -176,
                                     -176, -176, -176, -176, -176, -176, -176, -176, -175, -175, -175, -175,
                                     -175, -175, -175, -175, -175, -175, -175, -175, -174, -174, -174, -174,
                                     -174, -174, -174, -174, -174, -174, -174, -174, -173, -173, -173, -173,
                                     -173, -173, -173, -173, -173, -172, -172, -172, -172, -172, -172, -172,
                                     -172, -172, -172, -172, -172, -172, -172, -171, -171, -171, -171, -171,
                                     -171, -171, -171, -171, -170, -170, -170, -170, -170, -170, -170, -170,
                                     -170, -170, -170, -169, -169, -169, -169, -169, -169, -169, -169, -169,
                                     -169, -169, -169, -169, -169, -169, -169, -169, -169, -168, -168, -168,
                                     -168, -168, -168, -168, -168, -168, -168, -168, -167, -167, -167, -167,
                                     -167, -167, -167, -167, -166, -166, -166, -166, -166, -166, -166, -166,
                                     -166, -166, -166, -166, -166, -166, -166, -166, -165, -165, -165, -165,
                                     -165, -165, -165, -165, -165, -164, -164, -164, -164, -164, -164, -164,
                                     -164, -164, -163, -163, -163, -163, -163, -163, -163, -163, -162, -162,
                                     -162, -162, -162, -162, -162, -162, -162, -161, -161, -161, -161, -161,
                                     -161, -161, -161, -161, -161, -161, -161, -161, -160, -160, -160, -160,
                                     -160, -160, -160, -160, -160, -160, -160, -159, -159, -159, -159, -159,
                                     -159, -159, -159, -159, -159, -159, -159, -159, -159, -158, -158, -158,
                                     -158, -158, -158, -158, -158, -158, -158, -158, -158, -157, -157, -157,
                                     -157, -157, -157, -157, -157, -157, -156, -156, -156, -156, -156, -156,
                                     -156, -156, -156, -156, -156, -156, -156, -156, -155, -155, -155, -155,
                                     -155, -155, -155, -155, -155, -155, -155, -155, -155, -154, -154, -154,
                                     -154, -154, -154, -153, -153, -153, -153, -153, -153, -153, -153, -153,
                                     -153, -153, -153, -153, -152, -152, -152, -152, -152, -152, -152, -152,
                                     -152, -152, -152, -152, -152, -152, -152, -151, -151, -151, -151, -151,
                                     -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150,
                                     -150, -150, -150, -150, -149, -149, -149, -149, -149, -149, -149, -148,
                                     -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148,
                                     -147, -147, -147, -147, -147, -147, -147, -147, -147, -146, -146, -146,
                                     -146, -146, -146, -146, -146, -146, -146, -146, -146, -145, -145, -145,
                                     -145, -145, -145, -145, -145, -145, -144, -144, -144, -144, -144, -144,
                                     -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -142, -142,
                                     -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142,
                                     -142, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141,
                                     -140, -140, -140, -140, -140, -140, -140, -140, -140, -140, -139, -139,
                                     -139, -139, -139, -139, -139, -139, -138, -138, -138, -138, -138, -138,
                                     -138, -138, -138, -138, -137, -137, -137, -137, -137, -137, -137, -137,
                                     -137, -137, -137, -137, -137, -137, -137, -137, -136, -136, -136, -136,
                                     -136, -136, -136, -136, -136, -135, -135, -135, -135, -135, -135, -135,
                                     -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -134, -134,
                                     -134, -134, -134, -134, -134, -134, -134, -134, -133, -133, -133, -133,
                                     -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133,
                                     -133, -132, -132, -132, -132, -132, -132, -132, -132, -132, -131, -131,
                                     -131, -131, -131, -131, -131, -130, -130, -130, -130, -130, -130, -130,
                                     -130, -130, -130, -130, -130, -130, -129, -129, -129, -129, -129, -129,
                                     -129, -129, -129, -129, -129, -129, -129, -128, -128, -128, -128, -128,
                                     -128, -128, -128, -128, -128, -127, -127, -127, -127, -127, -127, -127,
                                     -127, -127, -127, -127, -126, -126, -126, -126, -126, -126, -125, -125,
                                     -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125,
                                     -125, -125, -125, -124, -124, -124, -124, -124, -124, -123, -123, -123,
                                     -123, -123, -123, -123, -123, -123, -122, -122, -122, -122, -122, -122,
                                     -122, -122, -122, -122, -122, -122, -121, -121, -121, -121, -121, -121,
                                     -121, -121, -121, -121, -121, -120, -120, -120, -120, -120, -120, -120,
                                     -120, -120, -120, -120, -120, -120, -120, -120, -120, -120, -120, -120,
                                     -120, -119, -119, -119, -119, -119, -119, -119, -119, -119, -119, -119,
                                     -119, -118, -118, -118, -118, -118, -118, -118, -118, -118, -118, -118,
                                     -117, -117, -117, -117, -117, -117, -117, -117, -117, -117, -117, -117,
                                     -116, -116, -116, -116, -116, -116, -116, -115, -115, -115, -115, -115,
                                     -115, -115, -115, -115, -115, -115, -115, -115, -114, -114, -114, -114,
                                     -114, -114, -114, -114, -114, -114, -114, -114, -114, -114, -114, -114,
                                     -113, -113, -113, -113, -113, -113, -113, -113, -112, -112, -112, -112,
                                     -112, -112, -112, -112, -112, -112, -112, -112, -112, -111, -111, -111,
                                     -111, -111, -111, -111, -111, -111, -111, -110, -110, -110, -110, -110,
                                     -110, -110, -110, -110, -110, -110, -110, -109, -109, -109, -109, -109,
                                     -109, -109, -109, -109, -109, -109, -109, -108, -108, -108, -108, -108,
                                     -108, -108, -108, -108, -108, -108, -107, -107, -107, -107, -107, -107,
                                     -107, -107, -107, -107, -106, -106, -106, -106, -106, -106, -106, -106,
                                     -106, -106, -106, -106, -105, -105, -105, -105, -105, -105, -105, -105,
                                     -105, -105, -105, -104, -104, -104, -104, -104, -104, -103, -103, -103,
                                     -103, -103, -103, -103, -103, -103, -103, -103, -103, -103, -103, -103,
                                     -103, -103, -102, -102, -102, -102, -102, -102, -102, -101, -101, -101,
                                     -101, -101, -101, -101, -101, -101, -101, -101, -101, -101, -101, -100,
                                     -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100,
                                     -100,  -99,  -99,  -99,  -99,  -99,  -99,  -98,  -98,  -98,  -98,  -98,
                                         -98,  -98,  -98,  -98,  -98,  -98,  -98,  -98,  -98,  -98,  -98,  -98,
                                         -98,  -97,  -97,  -97,  -97,  -97,  -97,  -97,  -97,  -97,  -97,  -97,
                                         -97,  -96,  -96,  -96,  -96,  -96,  -96,  -96,  -96,  -96,  -96,  -96,
                                         -95,  -95,  -95,  -95,  -95,  -95,  -95,  -95,  -95,  -95,  -95,  -95,
                                         -95,  -94,  -94,  -94,  -94,  -94,  -94,  -94,  -94,  -94,  -94,  -94,
                                         -94,  -94,  -94,  -93,  -93,  -93,  -93,  -93,  -93,  -93,  -93,  -93,
                                         -92,  -92,  -92,  -92,  -92,  -92,  -92,  -92,  -92,  -91,  -91,  -91,
                                         -91,  -91,  -91,  -91,  -91,  -91,  -91,  -91,  -91,  -91,  -91,  -90,
                                         -90,  -90,  -90,  -90,  -90,  -90,  -89,  -89,  -89,  -89,  -89,  -89,
                                         -89,  -89,  -89,  -89,  -89,  -89,  -89,  -89,  -88,  -88,  -88,  -88,
                                         -88,  -88,  -88,  -88,  -88,  -88,  -88,  -88,  -88,  -88,  -88,  -87,
                                         -87,  -87,  -87,  -87,  -87,  -87,  -87,  -87,  -87,  -87,  -87,  -86,
                                         -86,  -86,  -86,  -86,  -86,  -86,  -86,  -86,  -85,  -85,  -85,  -85,
                                         -85,  -85,  -85,  -85,  -85,  -85,  -85,  -84,  -84,  -84,  -84,  -84,
                                         -84,  -84,  -84,  -84,  -84,  -84,  -84,  -83,  -83,  -83,  -83,  -83,
                                         -83,  -82,  -82,  -82,  -82,  -82,  -82,  -82,  -82,  -82,  -82,  -82,
                                         -82,  -82,  -82,  -82,  -81,  -81,  -81,  -81,  -81,  -81,  -81,  -81,
                                         -81,  -81,  -81,  -81,  -81,  -81,  -80,  -80,  -80,  -80,  -80,  -80,
                                         -80,  -80,  -80,  -79,  -79,  -79,  -79,  -79,  -79,  -79,  -79,  -79,
                                         -79,  -78,  -78,  -78,  -78,  -78,  -78,  -78,  -78,  -78,  -78,  -77,
                                         -77,  -77,  -77,  -77,  -77,  -77,  -77,  -77,  -77,  -77,  -77,  -77,
                                         -77,  -77,  -76,  -76,  -76,  -76,  -76,  -76,  -76,  -76,  -76,  -76,
                                         -76,  -76,  -75,  -75,  -75,  -75,  -75,  -75,  -75,  -75,  -75,  -75,
                                         -75,  -75,  -75,  -75,  -75,  -75,  -75,  -74,  -74,  -74,  -74,  -74,
                                         -74,  -74,  -74,  -74,  -74,  -73,  -73,  -73,  -73,  -73,  -73,  -73,
                                         -73,  -73,  -73,  -72,  -72,  -72,  -72,  -72,  -72,  -72,  -71,  -71,
                                         -71,  -71,  -71,  -71,  -71,  -71,  -71,  -71,  -71,  -70,  -70,  -70,
                                         -70,  -70,  -70,  -70,  -70,  -70,  -70,  -70,  -70,  -70,  -70,  -70,
                                         -70,  -69,  -69,  -69,  -69,  -69,  -69,  -69,  -69,  -69,  -69,  -69,
                                         -69,  -69,  -68,  -68,  -68,  -68,  -68,  -68,  -68,  -68,  -68,  -68,
                                         -68,  -67,  -67,  -67,  -67,  -67,  -67,  -66,  -66,  -66,  -66,  -66,
                                         -66,  -66,  -66,  -66,  -66,  -66,  -66,  -66,  -66,  -66,  -66,  -65,
                                         -65,  -65,  -65,  -65,  -65,  -65,  -65,  -65,  -65,  -65,  -65,  -65,
                                         -65,  -64,  -64,  -64,  -64,  -64,  -64,  -64,  -64,  -64,  -64,  -64,
                                         -64,  -64,  -64,  -64,  -64,  -63,  -63,  -63,  -63,  -63,  -63,  -63,
                                         -62,  -62,  -62,  -62,  -62,  -62,  -62,  -61,  -61,  -61,  -61,  -61,
                                         -61,  -61,  -61,  -61,  -61,  -61,  -61,  -60,  -60,  -60,  -60,  -60,
                                         -60,  -60,  -60,  -60,  -60,  -60,  -60,  -60,  -60,  -60,  -60,  -59,
                                         -59,  -59,  -59,  -59,  -59,  -59,  -59,  -59,  -59,  -59,  -59,  -59,
                                         -59,  -59,  -59,  -59,  -59,  -59,  -59,  -59,  -58,  -58,  -58,  -58,
                                         -58,  -58,  -58,  -58,  -58,  -57,  -57,  -57,  -57,  -57,  -57,  -57,
                                         -57,  -57,  -57,  -57,  -56,  -56,  -56,  -56,  -56,  -56,  -56,  -56,
                                         -56,  -56,  -56,  -55,  -55,  -55,  -55,  -55,  -55,  -55,  -55,  -55,
                                         -55,  -55,  -55,  -55,  -55,  -54,  -54,  -54,  -54,  -54,  -54,  -54,
                                         -54,  -54,  -54,  -54,  -54,  -54,  -54,  -54,  -54,  -53,  -53,  -53,
                                         -53,  -52,  -52,  -52,  -52,  -52,  -52,  -52,  -52,  -52,  -52,  -52,
                                         -52,  -51,  -51,  -51,  -51,  -51,  -51,  -51,  -51,  -51,  -51,  -51,
                                         -51,  -51,  -51,  -50,  -50,  -50,  -50,  -50,  -50,  -50,  -50,  -50,
                                         -50,  -50,  -50,  -50,  -49,  -49,  -49,  -49,  -49,  -49,  -49,  -49,
                                         -49,  -49,  -49,  -49,  -49,  -49,  -49,  -49,  -49,  -48,  -48,  -48,
                                         -48,  -48,  -48,  -48,  -48,  -48,  -48,  -48,  -48,  -48,  -47,  -47,
                                         -47,  -47,  -47,  -47,  -47,  -47,  -47,  -47,  -46,  -46,  -46,  -46,
                                         -46,  -46,  -46,  -46,  -46,  -46,  -46,  -46,  -45,  -45,  -45,  -45,
                                         -45,  -45,  -45,  -45,  -45,  -45,  -44,  -44,  -44,  -44,  -44,  -44,
                                         -44,  -44,  -44,  -44,  -44,  -44,  -44,  -44,  -43,  -43,  -43,  -43,
                                         -43,  -43,  -43,  -43,  -43,  -43,  -43,  -42,  -42,  -42,  -42,  -42,
                                         -42,  -42,  -42,  -42,  -42,  -42,  -42,  -41,  -41,  -41,  -41,  -41,
                                         -41,  -41,  -41,  -41,  -41,  -41,  -41,  -40,  -40,  -40,  -40,  -40,
                                         -40,  -40,  -40,  -40,  -40,  -40,  -40,  -40,  -40,  -39,  -39,  -39,
                                         -39,  -39,  -39,  -39,  -39,  -39,  -39,  -39,  -39,  -39,  -39,  -39,
                                         -38,  -38,  -38,  -38,  -38,  -38,  -38,  -38,  -38,  -38,  -38,  -38,
                                         -37,  -37,  -37,  -37,  -37,  -37,  -37,  -37,  -37,  -37,  -37,  -37,
                                         -37,  -37,  -37,  -37,  -37,  -36,  -36,  -36,  -36,  -36,  -36,  -36,
                                         -36,  -36,  -36,  -35,  -35,  -35,  -35,  -35,  -35,  -35,  -35,  -35,
                                         -34,  -34,  -34,  -34,  -34,  -34,  -34,  -34,  -34,  -34,  -34,  -34,
                                         -34,  -34,  -34,  -34,  -34,  -34,  -33,  -33,  -33,  -33,  -33,  -33,
                                         -33,  -33,  -33,  -33,  -33,  -33,  -33,  -32,  -32,  -32,  -32,  -32,
                                         -32,  -32,  -32,  -32,  -32,  -32,  -32,  -32,  -32,  -32,  -31,  -31,
                                         -31,  -31,  -31,  -31,  -31,  -31,  -31,  -31,  -31,  -31,  -31,  -30,
                                         -30,  -30,  -30,  -30,  -30,  -30,  -30,  -30,  -30,  -30,  -30,  -30,
                                         -30,  -29,  -29,  -29,  -29,  -29,  -29,  -29,  -29,  -29,  -29,  -29,
                                         -29,  -29,  -29,  -28,  -28,  -28,  -28,  -28,  -28,  -28,  -28,  -27,
                                         -27,  -27,  -27,  -27,  -27,  -27,  -27,  -27,  -27,  -27,  -27,  -27,
                                         -27,  -27,  -27,  -26,  -26,  -26,  -26,  -26,  -26,  -25,  -25,  -25,
                                         -25,  -25,  -25,  -25,  -25,  -25,  -25,  -25,  -25,  -25,  -24,  -24,
                                         -24,  -24,  -24,  -24,  -23,  -23,  -23,  -23,  -23,  -23,  -23,  -23,
                                         -23,  -22,  -22,  -22,  -22,  -22,  -22,  -22,  -22,  -22,  -22,  -21,
                                         -21,  -21,  -21,  -21,  -21,  -21,  -21,  -21,  -21,  -21,  -21,  -21,
                                         -21,  -20,  -20,  -20,  -20,  -20,  -20,  -20,  -20,  -20,  -20,  -20,
                                         -20,  -20,  -20,  -20,  -19,  -19,  -19,  -19,  -19,  -19,  -19,  -19,
                                         -19,  -19,  -19,  -19,  -19,  -19,  -18,  -18,  -18,  -18,  -18,  -18,
                                         -18,  -18,  -18,  -18,  -18,  -18,  -18,  -18,  -18,  -18,  -18,  -17,
                                         -17,  -17,  -17,  -17,  -17,  -17,  -17,  -17,  -17,  -17,  -17,  -17,
                                         -17,  -17,  -17,  -17,  -17,  -17,  -17,  -17,  -16,  -16,  -16,  -16,
                                         -16,  -16,  -16,  -16,  -16,  -16,  -16,  -16,  -16,  -16,  -16,  -15,
                                         -15,  -15,  -15,  -15,  -15,  -15,  -15,  -15,  -15,  -15,  -15,  -15,
                                         -15,  -15,  -15,  -15,  -15,  -14,  -14,  -14,  -14,  -14,  -14,  -14,
                                         -14,  -14,  -14,  -13,  -13,  -13,  -13,  -13,  -13,  -13,  -13,  -13,
                                         -13,  -13,  -12,  -12,  -12,  -12,  -12,  -12,  -12,  -12,  -11,  -11,
                                         -11,  -11,  -11,  -11,  -11,  -11,  -11,  -11,  -11,  -10,  -10,  -10,
                                         -10,  -10,  -10,  -10,  -10,  -10,  -10,  -10,  -10,   -9,   -9,   -9,
                                             -9,   -9,   -9,   -9,   -9,   -9,   -9,   -9,   -9,   -8,   -8,   -8,
                                             -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -7,   -7,   -7,   -7,
                                             -7,   -7,   -7,   -7,   -7,   -6,   -6,   -6,   -6,   -6,   -5,   -5,
                                             -5,   -5,   -5,   -5,   -5,   -5,   -5,   -5,   -5,   -5,   -5,   -5,
                                             -4,   -4,   -4,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
                                             -3,   -3,   -3,   -3,   -3,   -2,   -2,   -2,   -2,   -2,   -1,   -1,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    2,
                                                 2,    2,    2,    2,    3,    3,    3,    3,    3,    3,    3,    3,
                                                 3,    3,    3,    3,    3,    3,    4,    4,    4,    5,    5,    5,
                                                 5,    5,    5,    5,    5,    5,    5,    5,    5,    5,    5,    6,
                                                 6,    6,    6,    6,    7,    7,    7,    7,    7,    7,    7,    7,
                                                 7,    8,    8,    8,    8,    8,    8,    8,    8,    8,    8,    8,
                                                 9,    9,    9,    9,    9,    9,    9,    9,    9,    9,    9,    9,
                                             10,   10,   10,   10,   10,   10,   10,   10,   10,   10,   10,   10,
                                             11,   11,   11,   11,   11,   11,   11,   11,   11,   11,   11,   12,
                                             12,   12,   12,   12,   12,   12,   12,   13,   13,   13,   13,   13,
                                             13,   13,   13,   13,   13,   13,   14,   14,   14,   14,   14,   14,
                                             14,   14,   14,   14,   15,   15,   15,   15,   15,   15,   15,   15,
                                             15,   15,   15,   15,   15,   15,   15,   15,   15,   15,   16,   16,
                                             16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,   16,
                                             16,   17,   17,   17,   17,   17,   17,   17,   17,   17,   17,   17,
                                             17,   17,   17,   17,   17,   17,   17,   17,   17,   17,   18,   18,
                                             18,   18,   18,   18,   18,   18,   18,   18,   18,   18,   18,   18,
                                             18,   18,   18,   19,   19,   19,   19,   19,   19,   19,   19,   19,
                                             19,   19,   19,   19,   19,   20,   20,   20,   20,   20,   20,   20,
                                             20,   20,   20,   20,   20,   20,   20,   20,   21,   21,   21,   21,
                                             21,   21,   21,   21,   21,   21,   21,   21,   21,   21,   22,   22,
                                             22,   22,   22,   22,   22,   22,   22,   22,   23,   23,   23,   23,
                                             23,   23,   23,   23,   23,   24,   24,   24,   24,   24,   24,   25,
                                             25,   25,   25,   25,   25,   25,   25,   25,   25,   25,   25,   25,
                                             26,   26,   26,   26,   26,   26,   27,   27,   27,   27,   27,   27,
                                             27,   27,   27,   27,   27,   27,   27,   27,   27,   27,   28,   28,
                                             28,   28,   28,   28,   28,   28,   29,   29,   29,   29,   29,   29,
                                             29,   29,   29,   29,   29,   29,   29,   29,   30,   30,   30,   30,
                                             30,   30,   30,   30,   30,   30,   30,   30,   30,   30,   31,   31,
                                             31,   31,   31,   31,   31,   31,   31,   31,   31,   31,   31,   32,
                                             32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,   32,
                                             32,   32,   33,   33,   33,   33,   33,   33,   33,   33,   33,   33,
                                             33,   33,   33,   34,   34,   34,   34,   34,   34,   34,   34,   34,
                                             34,   34,   34,   34,   34,   34,   34,   34,   34,   35,   35,   35,
                                             35,   35,   35,   35,   35,   35,   36,   36,   36,   36,   36,   36,
                                             36,   36,   36,   36,   37,   37,   37,   37,   37,   37,   37,   37,
                                             37,   37,   37,   37,   37,   37,   37,   37,   37,   38,   38,   38,
                                             38,   38,   38,   38,   38,   38,   38,   38,   38,   39,   39,   39,
                                             39,   39,   39,   39,   39,   39,   39,   39,   39,   39,   39,   39,
                                             40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,   40,
                                             40,   40,   41,   41,   41,   41,   41,   41,   41,   41,   41,   41,
                                             41,   41,   42,   42,   42,   42,   42,   42,   42,   42,   42,   42,
                                             42,   42,   43,   43,   43,   43,   43,   43,   43,   43,   43,   43,
                                             43,   44,   44,   44,   44,   44,   44,   44,   44,   44,   44,   44,
                                             44,   44,   44,   45,   45,   45,   45,   45,   45,   45,   45,   45,
                                             45,   46,   46,   46,   46,   46,   46,   46,   46,   46,   46,   46,
                                             46,   47,   47,   47,   47,   47,   47,   47,   47,   47,   47,   48,
                                             48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,   48,
                                             49,   49,   49,   49,   49,   49,   49,   49,   49,   49,   49,   49,
                                             49,   49,   49,   49,   49,   50,   50,   50,   50,   50,   50,   50,
                                             50,   50,   50,   50,   50,   50,   51,   51,   51,   51,   51,   51,
                                             51,   51,   51,   51,   51,   51,   51,   51,   52,   52,   52,   52,
                                             52,   52,   52,   52,   52,   52,   52,   52,   53,   53,   53,   53,
                                             54,   54,   54,   54,   54,   54,   54,   54,   54,   54,   54,   54,
                                             54,   54,   54,   54,   55,   55,   55,   55,   55,   55,   55,   55,
                                             55,   55,   55,   55,   55,   55,   56,   56,   56,   56,   56,   56,
                                             56,   56,   56,   56,   56,   57,   57,   57,   57,   57,   57,   57,
                                             57,   57,   57,   57,   58,   58,   58,   58,   58,   58,   58,   58,
                                             58,   59,   59,   59,   59,   59,   59,   59,   59,   59,   59,   59,
                                             59,   59,   59,   59,   59,   59,   59,   59,   59,   59,   60,   60,
                                             60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,   60,
                                             60,   60,   61,   61,   61,   61,   61,   61,   61,   61,   61,   61,
                                             61,   61,   62,   62,   62,   62,   62,   62,   62,   63,   63,   63,
                                             63,   63,   63,   63,   64,   64,   64,   64,   64,   64,   64,   64,
                                             64,   64,   64,   64,   64,   64,   64,   64,   65,   65,   65,   65,
                                             65,   65,   65,   65,   65,   65,   65,   65,   65,   65,   66,   66,
                                             66,   66,   66,   66,   66,   66,   66,   66,   66,   66,   66,   66,
                                             66,   66,   67,   67,   67,   67,   67,   67,   68,   68,   68,   68,
                                             68,   68,   68,   68,   68,   68,   68,   69,   69,   69,   69,   69,
                                             69,   69,   69,   69,   69,   69,   69,   69,   70,   70,   70,   70,
                                             70,   70,   70,   70,   70,   70,   70,   70,   70,   70,   70,   70,
                                             71,   71,   71,   71,   71,   71,   71,   71,   71,   71,   71,   72,
                                             72,   72,   72,   72,   72,   72,   73,   73,   73,   73,   73,   73,
                                             73,   73,   73,   73,   74,   74,   74,   74,   74,   74,   74,   74,
                                             74,   74,   75,   75,   75,   75,   75,   75,   75,   75,   75,   75,
                                             75,   75,   75,   75,   75,   75,   75,   76,   76,   76,   76,   76,
                                             76,   76,   76,   76,   76,   76,   76,   77,   77,   77,   77,   77,
                                             77,   77,   77,   77,   77,   77,   77,   77,   77,   77,   78,   78,
                                             78,   78,   78,   78,   78,   78,   78,   78,   79,   79,   79,   79,
                                             79,   79,   79,   79,   79,   79,   80,   80,   80,   80,   80,   80,
                                             80,   80,   80,   81,   81,   81,   81,   81,   81,   81,   81,   81,
                                             81,   81,   81,   81,   81,   82,   82,   82,   82,   82,   82,   82,
                                             82,   82,   82,   82,   82,   82,   82,   82,   83,   83,   83,   83,
                                             83,   83,   84,   84,   84,   84,   84,   84,   84,   84,   84,   84,
                                             84,   84,   85,   85,   85,   85,   85,   85,   85,   85,   85,   85,
                                             85,   86,   86,   86,   86,   86,   86,   86,   86,   86,   87,   87,
                                             87,   87,   87,   87,   87,   87,   87,   87,   87,   87,   88,   88,
                                             88,   88,   88,   88,   88,   88,   88,   88,   88,   88,   88,   88,
                                             88,   89,   89,   89,   89,   89,   89,   89,   89,   89,   89,   89,
                                             89,   89,   89,   90,   90,   90,   90,   90,   90,   90,   91,   91,
                                             91,   91,   91,   91,   91,   91,   91,   91,   91,   91,   91,   91,
                                             92,   92,   92,   92,   92,   92,   92,   92,   92,   93,   93,   93,
                                             93,   93,   93,   93,   93,   93,   94,   94,   94,   94,   94,   94,
                                             94,   94,   94,   94,   94,   94,   94,   94,   95,   95,   95,   95,
                                             95,   95,   95,   95,   95,   95,   95,   95,   95,   96,   96,   96,
                                             96,   96,   96,   96,   96,   96,   96,   96,   97,   97,   97,   97,
                                             97,   97,   97,   97,   97,   97,   97,   97,   98,   98,   98,   98,
                                             98,   98,   98,   98,   98,   98,   98,   98,   98,   98,   98,   98,
                                             98,   98,   99,   99,   99,   99,   99,   99,  100,  100,  100,  100,
                                         100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  101,  101,
                                         101,  101,  101,  101,  101,  101,  101,  101,  101,  101,  101,  101,
                                         102,  102,  102,  102,  102,  102,  102,  103,  103,  103,  103,  103,
                                         103,  103,  103,  103,  103,  103,  103,  103,  103,  103,  103,  103,
                                         104,  104,  104,  104,  104,  104,  105,  105,  105,  105,  105,  105,
                                         105,  105,  105,  105,  105,  106,  106,  106,  106,  106,  106,  106,
                                         106,  106,  106,  106,  106,  107,  107,  107,  107,  107,  107,  107,
                                         107,  107,  107,  108,  108,  108,  108,  108,  108,  108,  108,  108,
                                         108,  108,  109,  109,  109,  109,  109,  109,  109,  109,  109,  109,
                                         109,  109,  110,  110,  110,  110,  110,  110,  110,  110,  110,  110,
                                         110,  110,  111,  111,  111,  111,  111,  111,  111,  111,  111,  111,
                                         112,  112,  112,  112,  112,  112,  112,  112,  112,  112,  112,  112,
                                         112,  113,  113,  113,  113,  113,  113,  113,  113,  114,  114,  114,
                                         114,  114,  114,  114,  114,  114,  114,  114,  114,  114,  114,  114,
                                         114,  115,  115,  115,  115,  115,  115,  115,  115,  115,  115,  115,
                                         115,  115,  116,  116,  116,  116,  116,  116,  116,  117,  117,  117,
                                         117,  117,  117,  117,  117,  117,  117,  117,  117,  118,  118,  118,
                                         118,  118,  118,  118,  118,  118,  118,  118,  119,  119,  119,  119,
                                         119,  119,  119,  119,  119,  119,  119,  119,  120,  120,  120,  120,
                                         120,  120,  120,  120,  120,  120,  120,  120,  120,  120,  120,  120,
                                         120,  120,  120,  120,  121,  121,  121,  121,  121,  121,  121,  121,
                                         121,  121,  121,  122,  122,  122,  122,  122,  122,  122,  122,  122,
                                         122,  122,  122,  123,  123,  123,  123,  123,  123,  123,  123,  123,
                                         124,  124,  124,  124,  124,  124,  125,  125,  125,  125,  125,  125,
                                         125,  125,  125,  125,  125,  125,  125,  125,  125,  125,  125,  126,
                                         126,  126,  126,  126,  126,  127,  127,  127,  127,  127,  127,  127,
                                         127,  127,  127,  127,  128,  128,  128,  128,  128,  128,  128,  128,
                                         128,  128,  129,  129,  129,  129,  129,  129,  129,  129,  129,  129,
                                         129,  129,  129,  130,  130,  130,  130,  130,  130,  130,  130,  130,
                                         130,  130,  130,  130,  131,  131,  131,  131,  131,  131,  131,  132,
                                         132,  132,  132,  132,  132,  132,  132,  132,  133,  133,  133,  133,
                                         133,  133,  133,  133,  133,  133,  133,  133,  133,  133,  133,  133,
                                         133,  134,  134,  134,  134,  134,  134,  134,  134,  134,  134,  135,
                                         135,  135,  135,  135,  135,  135,  135,  135,  135,  135,  135,  135,
                                         135,  135,  135,  135,  136,  136,  136,  136,  136,  136,  136,  136,
                                         136,  137,  137,  137,  137,  137,  137,  137,  137,  137,  137,  137,
                                         137,  137,  137,  137,  137,  138,  138,  138,  138,  138,  138,  138,
                                         138,  138,  138,  139,  139,  139,  139,  139,  139,  139,  139,  140,
                                         140,  140,  140,  140,  140,  140,  140,  140,  140,  141,  141,  141,
                                         141,  141,  141,  141,  141,  141,  141,  141,  142,  142,  142,  142,
                                         142,  142,  142,  142,  142,  142,  142,  142,  142,  142,  142,  143,
                                         143,  143,  143,  143,  143,  143,  143,  143,  143,  144,  144,  144,
                                         144,  144,  144,  145,  145,  145,  145,  145,  145,  145,  145,  145,
                                         146,  146,  146,  146,  146,  146,  146,  146,  146,  146,  146,  146,
                                         147,  147,  147,  147,  147,  147,  147,  147,  147,  148,  148,  148,
                                         148,  148,  148,  148,  148,  148,  148,  148,  148,  148,  149,  149,
                                         149,  149,  149,  149,  149,  150,  150,  150,  150,  150,  150,  150,
                                         150,  150,  150,  150,  150,  150,  150,  150,  150,  151,  151,  151,
                                         151,  151,  152,  152,  152,  152,  152,  152,  152,  152,  152,  152,
                                         152,  152,  152,  152,  152,  153,  153,  153,  153,  153,  153,  153,
                                         153,  153,  153,  153,  153,  153,  154,  154,  154,  154,  154,  154,
                                         155,  155,  155,  155,  155,  155,  155,  155,  155,  155,  155,  155,
                                         155,  156,  156,  156,  156,  156,  156,  156,  156,  156,  156,  156,
                                         156,  156,  156,  157,  157,  157,  157,  157,  157,  157,  157,  157,
                                         158,  158,  158,  158,  158,  158,  158,  158,  158,  158,  158,  158,
                                         159,  159,  159,  159,  159,  159,  159,  159,  159,  159,  159,  159,
                                         159,  159,  160,  160,  160,  160,  160,  160,  160,  160,  160,  160,
                                         160,  161,  161,  161,  161,  161,  161,  161,  161,  161,  161,  161,
                                         161,  161,  162,  162,  162,  162,  162,  162,  162,  162,  162,  163,
                                         163,  163,  163,  163,  163,  163,  163,  164,  164,  164,  164,  164,
                                         164,  164,  164,  164,  165,  165,  165,  165,  165,  165,  165,  165,
                                         165,  166,  166,  166,  166,  166,  166,  166,  166,  166,  166,  166,
                                         166,  166,  166,  166,  166,  167,  167,  167,  167,  167,  167,  167,
                                         167,  168,  168,  168,  168,  168,  168,  168,  168,  168,  168,  168,
                                         169,  169,  169,  169,  169,  169,  169,  169,  169,  169,  169,  169,
                                         169,  169,  169,  169,  169,  169,  170,  170,  170,  170,  170,  170,
                                         170,  170,  170,  170,  170,  171,  171,  171,  171,  171,  171,  171,
                                         171,  171,  172,  172,  172,  172,  172,  172,  172,  172,  172,  172,
                                         172,  172,  172,  172,  173,  173,  173,  173,  173,  173,  173,  173,
                                         173,  174,  174,  174,  174,  174,  174,  174,  174,  174,  174,  174,
                                         174,  175,  175,  175,  175,  175,  175,  175,  175,  175,  175,  175,
                                         175,  176,  176,  176,  176,  176,  176,  176,  176,  176,  176,  176,
                                         176,  177,  177,  177,  177,  177,  177,  177,  177,  177,  177,  177,
                                         177,  177,  177,  177,  177,  178,  178,  178,  178,  178,  178,  178,
                                         179,  179,  179,  179,  179,  179,  179,  179,  179,  180,  180,  180,
                                         180,  180,  180,  180,  180,  180,  180,  180,  181,  181,  181,  181,
                                         181,  181,  181,  181,  182,  182,  182,  182,  182,  182,  182,  182,
                                         182,  182,  183,  183,  183,  183,  183,  183,  183,  184,  184,  184,
                                         184,  184,  184,  184,  184,  184,  184,  184,  184,  184,  184,  184,
                                         184,  185,  185,  185,  185,  185,  185,  185,  185,  185,  185,  185,
                                         185,  186,  186,  186,  186,  186,  186,  186,  186,  186,  186,  186,
                                         186,  186,  187,  187,  187,  187,  187,  187,  187,  187,  187,  187,
                                         187,  187,  188,  188,  188,  188,  188,  188,  188,  188,  188,  188,
                                         188,  189,  189,  189,  189,  189,  189,  189,  189,  189,  189,  190,
                                         190,  190,  190,  190,  190,  190,  190,  190,  190,  191,  191,  191,
                                         191,  191,  191,  191,  191,  191,  192,  192,  192,  192,  192,  192,
                                         192,  192,  192,  193,  193,  193,  193,  193,  193,  193,  193,  194,
                                         194,  194,  194,  194,  194,  194,  194,  194,  194,  194,  195,  195,
                                         195,  195,  195,  195,  195,  195,  195,  195,  195,  195,  196,  196,
                                         196,  196,  196,  196,  197,  197,  197,  197,  197,  197,  197,  197,
                                         197,  197,  197,  197,  197,  197,  198,  198,  198,  198,  198,  198,
                                         198,  199,  199,  199,  199,  199,  199,  199,  199,  199,  199,  199,
                                         199,  199,  199,  199,  199,  200,  200,  200,  200,  200,  200,  200,
                                         200,  200,  200,  201,  201,  201,  201,  201,  201,  201,  201,  201,
                                         201,  201,  202,  202,  202,  202,  202,  202,  202,  202,  202,  202,
                                         203,  203,  203,  203,  203,  203,  203,  203,  203,  203,  203,  203,
                                         204,  204,  204,  204,  204,  204,  204,  204,  204,  204,  204,  205,
                                         205,  205,  205,  205,  205,  205,  205,  205,  206,  206,  206,  206,
                                         206,  206,  206,  206,  206,  206,  206,  206,  207,  207,  207,  207,
                                         207,  207,  207,  207,  207,  207,  207,  207,  208,  208,  208,  208,
                                         208,  208,  208,  208,  208,  208,  209,  209,  209,  209,  209,  209,
                                         209,  209,  210,  210,  210,  210,  210,  210,  210,  210,  210,  211,
                                         211,  211,  211,  211,  211,  211,  211,  212,  212,  212,  212,  212,
                                         212,  212,  212,  212,  212,  212,  212,  213,  213,  213,  213,  213,
                                         213,  213,  213,  213,  213,  213,  214,  214,  214,  214,  214,  214,
                                         214,  214,  214,  215,  215,  215,  215,  215,  215,  215,  215,  215,
                                         215,  215,  215,  216,  216,  216,  216,  216,  216,  216,  216,  216,
                                         216,  216,  216,  216,  217,  217,  217,  217,  217,  217,  217,  217,
                                         217,  218,  218,  218,  218,  218,  218,  218,  218,  218,  218,  218,
                                         218,  218,  219,  219,  219,  219,  219,  219,  219,  220,  220,  220,
                                         220,  220,  220,  220,  220,  220,  221,  221,  221,  221,  221,  221,
                                         221,  221,  221,  221,  221,  221,  221,  222,  222,  222,  222,  222,
                                         223,  223,  223,  223,  223,  223,  223,  223,  223,  223,  223,  224,
                                         224,  224,  224,  224,  224,  224,  224,  224,  225,  225,  225,  225,
                                         225,  225,  225,  225,  226,  226,  226,  226,  226,  226,  226,  226,
                                         226,  226,  226,  226,  226,  226,  227,  227,  227,  227,  227,  227,
                                         228,  228,  228,  228,  228,  228,  228,  228,  228,  228,  228,  229,
                                         229,  229,  229,  229,  229,  229,  229,  229,  230,  230,  230,  230,
                                         230,  230,  230,  230,  230,  230,  230,  230,  230,  230,  230,  231,
                                         231,  231,  231,  231,  231,  231,  231,  231,  231,  231,  232,  232,
                                         232,  232,  232,  232,  232,  232,  232,  233,  233,  233,  233,  233,
                                         233,  233,  233,  233,  233,  233,  234,  234,  234,  234,  234,  234,
                                         234,  234,  234,  234,  234,  234,  234,  234,  234,  234,  235,  235,
                                         235,  235,  235,  235,  235,  235,  235,  235,  236,  236,  236,  236,
                                         236,  236,  236,  236,  236,  236,  236,  236,  237,  237,  237,  237,
                                         237,  237,  237,  237,  238,  238,  238,  238,  238,  238,  238,  239,
                                         239,  239,  239,  239,  239,  239,  239,  239,  240,  240,  240,  240,
                                         240,  240,  240,  241,  241,  241,  241,  241,  241,  241,  241,  241,
                                         241,  241,  242,  242,  242,  242,  242,  242,  242,  242,  243,  243,
                                         243,  243,  243,  243,  243,  243,  244,  244,  244,  244,  244,  244,
                                         244,  244,  244,  244,  245,  245,  245,  245,  245,  245,  245,  245,
                                         245,  246,  246,  246,  246,  246,  246,  246,  246,  246,  246,  246,
                                         247,  247,  247,  247,  247,  247,  247,  247,  247,  247,  247,  247,
                                         247,  248,  248,  248,  248,  248,  248,  248,  248,  249,  249,  249,
                                         249,  249,  249,  250,  250,  250,  250,  250,  251,  251,  251,  251,
                                         251,  251,  251,  251,  251,  251,  251,  251,  251,  251,  251,  252,
                                         252,  252,  252,  252,  252,  252,  252,  253,  253,  253,  253,  253,
                                         253,  253,  253,  253,  253,  253,  254,  254,  254,  254,  254,  254,
                                         254,  254,  255,  255,  255,  255,  255,  255,  255,  255,  255,  255,
                                         255,  255,  255,  256,  256,  256,  256,  256,  256,  256,  256,  256,
                                         256,  256,  256,  256,  257,  257,  257,  257,  257,  257,  257,  257,
                                         257,  257,  257,  258,  258,  258,  258,  258,  258,  258,  258,  258,
                                         258,  258,  259,  259,  259,  259,  259,  260,  260,  260,  260,  260,
                                         260,  260,  260,  260,  260,  261,  261,  261,  261,  261,  261,  261,
                                         261,  261,  261,  261,  262,  262,  262,  262,  262,  262,  262,  262,
                                         263,  263,  263,  263,  263,  263,  263,  263,  263,  263,  263,  264,
                                         264,  264,  264,  264,  264,  264,  264,  264,  265,  265,  265,  265,
                                         265,  265,  266,  266,  266,  266,  266,  266,  266,  266,  266,  266,
                                         266,  266,  266,  266,  267,  267,  267,  267,  267,  267,  267,  267,
                                         268,  268,  268,  268,  268,  268,  268,  268,  268,  268,  269,  269,
                                         269,  269,  269,  269,  269,  269,  270,  270,  270,  270,  270,  270,
                                         270,  270,  270,  271,  271,  271,  271,  271,  271,  271,  271,  271,
                                         271,  271,  272,  272,  272,  272,  272,  272,  272,  272,  272,  272,
                                         272,  272,  272,  273,  273,  273,  273,  273,  273,  273,  273,  273,
                                         273,  273,  273,  273,  274,  274,  274,  274,  274,  274,  274,  274,
                                         275,  275,  275,  275,  275,  275,  275,  275,  275,  275,  275,  275,
                                         275,  275,  276,  276,  276,  276,  276,  276,  276,  276,  276,  276,
                                         276,  277,  277,  277,  277,  277,  277,  277,  277,  277,  278,  278,
                                         278,  278,  278,  278,  278,  278,  278,  278,  278,  278,  278,  278,
                                         279,  279,  279,  279,  279,  279,  280,  280,  280,  280,  280,  280,
                                         280,  280,  280,  281,  281,  281,  281,  281,  281,  281,  282,  282,
                                         282,  282,  282,  282,  283,  283,  283,  283,  283,  283,  283,  283,
                                         283,  284,  284,  284,  284,  284,  284,  285,  285,  285,  285,  285,
                                         285,  285,  285,  285,  285,  285,  285,  285,  285,  285,  286,  286,
                                         286,  286,  286,  286,  286,  286,  286,  286,  287,  287,  287,  287,
                                         287,  287,  288,  288,  288,  288,  288,  288,  288,  288,  288,  288,
                                         288,  288,  288,  289,  289,  289,  289,  289,  289,  289,  289,  289,
                                         290,  290,  290,  290,  290,  291,  291,  291,  291,  291,  291,  291,
                                         291,  291,  291,  291,  291,  291,  291,  291,  291,  291,  291,  291,
                                         292,  292,  292,  292,  292,  292,  292,  292,  293,  293,  293,  293,
                                         293,  293,  293,  293,  294,  294,  294,  294,  294,  294,  294,  294,
                                         294,  294,  294,  294,  294,  295,  295,  295,  295,  295,  295,  295,
                                         295,  295,  295,  295,  295,  296,  296,  296,  296,  296,  296,  296,
                                         296,  296,  296,  296,  297,  297,  297,  297,  297,  297,  297,  297,
                                         297,  297,  298,  298,  298,  298,  298,  298,  298,  298,  298,  299,
                                         299,  299,  299,  299,  299,  300,  300,  300,  300,  300,  300,  300,
                                         300,  301,  301,  301,  301,  301,  301,  301,  301,  302,  302,  302,
                                         302,  302,  302,  302,  302,  302,  302,  303,  303,  303,  303,  303,
                                         303,  303,  303,  303,  303,  303,  304,  304,  304,  304,  304,  305,
                                         305,  305,  305,  305,  305,  305,  305,  305,  305,  305,  305,  305,
                                         305,  305,  306,  306,  306,  306,  306,  306,  306,  306,  307,  307,
                                         307,  307,  307,  307,  307,  307,  307,  307,  307,  307,  307,  308,
                                         308,  308,  308,  308,  308,  308,  308,  308,  308,  308,  309,  309,
                                         309,  309,  309,  309,  309,  309,  309,  310,  310,  310,  310,  310,
                                         310,  310,  310,  310,  310,  311,  311,  311,  311,  311,  312,  312,
                                         312,  312,  312,  312,  312,  312,  312,  313,  313,  313,  313,  313,
                                         313,  313,  314,  314,  314,  314,  314,  314,  314,  314,  314,  315,
                                         315,  315,  315,  315,  315,  315,  315,  315,  315,  315,  315,  316,
                                         316,  316,  316,  316,  316,  316,  316,  316,  316,  316,  316,  316,
                                         316,  316,  317,  317,  317,  317,  317,  317,  317,  317,  318,  318,
                                         318,  318,  318,  318,  318,  318,  318,  318,  319,  319,  319,  319,
                                         319,  319,  319,  319,  319,  319,  319,  320,  320,  320,  320,  320,
                                         320,  321,  321,  321,  321,  321,  321,  321,  321,  321,  321,  321,
                                         321,  322,  322,  322,  322,  322,  322,  322,  322,  322,  323,  323,
                                         323,  323,  323,  323,  323,  323,  323,  323,  324,  324,  324,  324,
                                         324,  324,  324,  324,  324,  325,  325,  325,  325,  325,  325,  325,
                                         325,  326,  326,  326,  326,  326,  326,  326,  326,  326,  326,  326,
                                         327,  327,  327,  327,  327,  327,  327,  327,  328,  328,  328,  328,
                                         328,  328,  328,  328,  328,  328,  328,  328,  328,  328,  329,  329,
                                         329,  329,  329,  329,  329,  330,  330,  330,  330,  330,  330,  330,
                                         331,  331,  331,  331,  331,  331,  331,  331,  332,  332,  332,  332,
                                         332,  332,  332,  332,  332,  332,  332,  332,  332,  333,  333,  333,
                                         333,  333,  333,  333,  333,  333,  333,  333,  333,  334,  334,  334,
                                         334,  334,  334,  334,  334,  334,  334,  334,  334,  335,  335,  335,
                                         335,  335,  335,  335,  335,  335,  335,  336,  336,  336,  336,  336,
                                         336,  336,  336,  336,  336,  336,  337,  337,  337,  337,  337,  337,
                                         337,  337,  337,  337,  338,  338,  338,  338,  338,  338,  338,  338,
                                         338,  338,  339,  339,  339,  339,  339,  339,  339,  339,  340,  340,
                                         340,  340,  340,  340,  340,  340,  340,  340,  340,  340,  341,  341,
                                         341,  341,  341,  341,  341,  341,  341,  341,  341,  342,  342,  342,
                                         342,  342,  343,  343,  343,  343,  343,  343,  343,  344,  344,  344,
                                         344,  344,  344,  344,  344,  345,  345,  345,  345,  345,  345,  345,
                                         345,  346,  346,  346,  346,  346,  346,  346,  346,  346,  347,  347,
                                         347,  347,  347,  347,  347,  347,  347,  347,  348,  348,  348,  348,
                                         348,  348,  348,  349,  349,  349,  349,  349,  349,  349,  349,  349,
                                         349,  349,  349,  350,  350,  350,  350,  350,  350,  350,  350,  350,
                                         351,  351,  351,  351,  351,  351,  351,  351,  351,  351,  351,  352,
                                         352,  352,  352,  352,  352,  352,  352,  353,  353,  353,  353,  353,
                                         353,  353,  353,  353,  354,  354,  354,  354,  354,  354,  354,  354,
                                         354,  354,  355,  355,  355,  355,  355,  355,  355,  355,  355,  355,
                                         355,  355,  355,  355,  355,  356,  356,  356,  356,  356,  356,  356,
                                         356,  356,  356,  356,  357,  357,  357,  357,  357,  357,  358,  358,
                                         358,  358,  358,  358,  358,  358,  358,  358,  358,  359,  359,  359,
                                         359,  359,  359,  360,  360,  360,  360,  360,  361,  361,  361,  361,
                                         361,  361,  361,  361,  361,  361,  361,  361,  361,  361,  362,  362,
                                         362,  362,  362,  363,  363,  363,  363,  363,  363,  363,  363,  363,
                                         363,  364,  364,  364,  364,  364,  364,  364,  365,  365,  365,  365,
                                         365,  365,  366,  366,  366,  366,  366,  366,  366,  366,  366,  366,
                                         366,  366,  367,  367,  367,  367,  367,  367,  367,  367,  367,  367,
                                         367,  367,  367,  368,  368,  368,  368,  368,  368,  368,  368,  368,
                                         368,  368,  368,  368,  369,  369,  369,  369,  369,  369,  369,  369,
                                         369,  369,  369,  369,  369,  369,  370,  370,  370,  370,  370,  370,
                                         370,  370,  371,  371,  371,  371,  371,  371,  371,  371,  371,  371,
                                         372,  372,  372,  372,  372,  372,  372,  372,  372,  372,  372,  372,
                                         373,  373,  373,  373,  373,  373,  373,  373,  373,  374,  374,  374,
                                         374,  374,  374,  374,  374,  375,  375,  375,  375,  375,  375,  375,
                                         376,  376,  376,  376,  376,  376,  376,  376,  376,  376,  376,  376,
                                         377,  377,  377,  377,  377,  378,  378,  378,  378,  378,  378,  378,
                                         378,  378,  378,  378,  378,  378,  379,  379,  379,  379,  379,  380,
                                         380,  380,  380,  380,  380,  381,  381,  381,  381,  381,  382,  382,
                                         382,  382,  382,  382,  382,  382,  382,  382,  382,  383,  383,  383,
                                         383,  383,  383,  383,  383,  383,  383,  384,  384,  384,  384,  384,
                                         384,  384,  384,  384,  384,  385,  385,  385,  385,  385,  385,  385,
                                         385,  385,  385,  385,  386,  386,  386,  386,  386,  386,  386,  386,
                                         386,  386,  386,  386,  387,  387,  387,  387,  387,  387,  387,  387,
                                         387,  387,  388,  388,  388,  388,  388,  388,  388,  388,  388,  388,
                                         389,  389,  389,  389,  389,  389,  389,  389,  389,  389,  389,  389,
                                         390,  390,  390,  390,  390,  390,  391,  391,  391,  391,  391,  391,
                                         391,  391,  391,  392,  392,  392,  392,  392,  392,  392,  392,  392,
                                         393,  393,  393,  393,  393,  393,  393,  393,  394,  394,  394,  394,
                                         394,  394,  394,  394,  394,  394,  394,  394,  394,  394,  395,  395,
                                         395,  395,  395,  396,  396,  396,  396,  396,  396,  396,  396,  397,
                                         397,  397,  397,  397,  397,  397,  397,  397,  397,  398,  398,  398,
                                         398,  398,  398,  398,  398,  398,  398,  399,  399,  399,  400,  400,
                                         400,  400,  400,  400,  401,  401,  401,  401,  401,  401,  401,  401,
                                         402,  402,  402,  402,  402,  402,  402,  403,  403,  403,  403,  403,
                                         403,  403,  403,  403,  403,  403,  403,  403,  403,  403,  403,  403,
                                         403,  403,  403,  403,  404,  404,  404,  404,  404,  404,  404,  404,
                                         405,  405,  405,  405,  405,  405,  405,  405,  405,  405,  405,  406,
                                         406,  406,  406,  406,  406,  406,  406,  406,  406,  406,  406,  407,
                                         407,  407,  407,  407,  407,  408,  408,  408,  408,  408,  408,  408,
                                         408,  409,  409,  409,  409,  409,  410,  410,  410,  410,  410,  410,
                                         410,  410,  410,  411,  411,  411,  411,  411,  411,  411,  412,  412,
                                         412,  412,  412,  412,  412,  412,  412,  412,  412,  412,  413,  413,
                                         413,  413,  413,  413,  413,  413,  413,  414,  414,  414,  414,  414,
                                         414,  414,  415,  415,  415,  415,  415,  415,  415,  415,  415,  415,
                                         416,  416,  416,  416,  416,  416,  416,  416,  416,  417,  417,  417,
                                         417,  417,  417,  417,  418,  418,  418,  418,  419,  419,  419,  419,
                                         419,  419,  419,  420,  420,  420,  420,  420,  420,  420,  420,  420,
                                         420,  420,  420,  421,  421,  421,  421,  421,  421,  421,  421,  421,
                                         421,  421,  422,  422,  422,  422,  422,  423,  423,  423,  423,  423,
                                         423,  423,  423,  424,  424,  424,  424,  424,  424,  424,  424,  425,
                                         425,  425,  425,  425,  425,  425,  425,  425,  425,  426,  426,  426,
                                         426,  426,  426,  427,  427,  427,  427,  427,  427,  427,  427,  427,
                                         427,  427,  428,  428,  428,  428,  428,  428,  428,  428,  428,  428,
                                         429,  429,  429,  429,  429,  429,  430,  430,  430,  430,  430,  430,
                                         430,  430,  430,  431,  431,  431,  431,  431,  431,  431,  431,  431,
                                         431,  431,  431,  432,  432,  432,  432,  432,  432,  432,  432,  432,
                                         433,  433,  433,  433,  433,  433,  433,  433,  433,  433,  433,  434,
                                         434,  434,  434,  434,  434,  435,  435,  435,  435,  435,  435,  435,
                                         435,  436,  436,  436,  436,  436,  436,  436,  436,  436,  436,  437,
                                         437,  437,  437,  437,  437,  437,  437,  438,  438,  438,  438,  438,
                                         438,  438,  438,  438,  438,  439,  439,  439,  440,  440,  440,  440,
                                         440,  440,  440,  441,  441,  441,  441,  441,  441,  441,  441,  441,
                                         441,  441,  442,  442,  442,  442,  442,  442,  442,  443,  443,  443,
                                         443,  443,  443,  443,  443,  443,  443,  443,  444,  444,  444,  444,
                                         444,  444,  444,  444,  445,  445,  445,  445,  445,  445,  445,  445,
                                         445,  445,  445,  446,  446,  446,  446,  446,  446,  446,  446,  446,
                                         446,  446,  446,  446,  447,  447,  447,  447,  447,  447,  447,  448,
                                         448,  448,  448,  448,  448,  449,  449,  449,  449,  449,  449,  449,
                                         449,  449,  449,  450,  450,  450,  450,  450,  450,  450,  450,  450,
                                         451,  451,  451,  451,  451,  451,  452,  452,  452,  452,  452,  452,
                                         452,  452,  452,  453,  453,  453,  453,  453,  453,  453,  453,  453,
                                         453,  454,  454,  454,  454,  454,  454,  454,  455,  455,  455,  455,
                                         455,  456,  456,  456,  456,  457,  457,  457,  457,  457,  457,  458,
                                         458,  458,  458,  458,  458,  458,  458,  458,  459,  459,  459,  459,
                                         459,  460,  460,  460,  460,  460,  460,  460,  460,  460,  461,  461,
                                         461,  461,  461,  461,  461,  461,  461,  461,  461,  461,  461,  461,
                                         461,  461,  462,  462,  462,  462,  462,  462,  462,  462,  462,  463,
                                         463,  463,  463,  463,  463,  463,  463,  463,  463,  464,  464,  464,
                                         464,  464,  464,  464,  464,  464,  464,  464,  464,  465,  465,  465,
                                         465,  465,  466,  466,  466,  466,  466,  466,  466,  466,  467,  467,
                                         467,  468,  468,  468,  468,  468,  468,  468,  468,  468,  468,  468,
                                         469,  469,  469,  469,  469,  469,  469,  470,  470,  470,  470,  470,
                                         471,  471,  471,  471,  471,  471,  471,  471,  471,  471,  471,  471,
                                         472,  472,  472,  472,  472,  472,  472,  473,  473,  473,  473,  473,
                                         473,  473,  474,  474,  474,  474,  475,  475,  475,  475,  475,  475,
                                         475,  475,  475,  475,  476,  476,  476,  476,  476,  476,  476,  476,
                                         476,  476,  476,  477,  477,  477,  477,  478,  478,  478,  478,  478,
                                         478,  478,  478,  478,  478,  478,  479,  479,  479,  479,  479,  479,
                                         479,  479,  479,  479,  479,  480,  480,  480,  480,  480,  480,  480,
                                         480,  481,  481,  481,  481,  481,  481,  481,  481,  481,  481,  482,
                                         482,  482,  482,  482,  482,  482,  482,  482,  483,  483,  483,  483,
                                         483,  483,  484,  484,  484,  484,  484,  484,  485,  485,  485,  485,
                                         485,  485,  485,  485,  485,  486,  486,  486,  486,  486,  486,  486,
                                         487,  487,  487,  487,  488,  488,  488,  488,  488,  488,  488,  488,
                                         488,  488,  488,  488,  488,  489,  489,  489,  489,  489,  489,  489,
                                         490,  490,  490,  490,  490,  490,  490,  490,  490,  490,  490,  491,
                                         491,  491,  491,  491,  491,  491,  491,  491,  491,  491,  491,  492,
                                         492,  492,  492,  492,  493,  493,  493,  493,  493,  493,  493,  493,
                                         494,  494,  494,  495,  495,  495,  495,  495,  495,  495,  495,  495,
                                         496,  496,  496,  496,  496,  496,  496,  496,  496,  497,  497,  497,
                                         497,  497,  497,  497,  497,  497,  497,  498,  498,  498,  498,  498,
                                         498,  498,  499,  499,  499,  499,  499,  499,  500,  500,  500,  500,
                                         500,  500,  500,  500,  500,  501,  501,  501,  501,  501,  501,  501,
                                         501,  502,  502,  502,  502,  502,  502,  502,  502,  502,  502,  502,
                                         503,  503,  503,  503,  503,  503,  503,  503,  503,  504,  504,  504,
                                         504,  504,  504,  504,  504,  504,  505,  505,  505,  505,  505,  505,
                                         505,  505,  505,  506,  506,  506,  506,  506,  506,  506,  506,  506,
                                         507,  507,  507,  507,  507,  508,  508,  508,  508,  508,  508,  508,
                                         508,  508,  509,  509,  509,  509,  509,  509,  509,  509,  509,  509,
                                         510,  510,  510,  510,  510,  510,  510,  511,  511,  511,  511,  511,
                                         511,  511,  511,  512,  512,  512,  512,  512,  512,  513,  513,  513,
                                         513,  513,  513,  513,  514,  514,  514,  514,  515,  515,  515,  515,
                                         515,  515,  515,  515,  515,  516,  516,  516,  517,  517,  517,  517,
                                         517,  517,  517,  517,  517,  517,  518,  518,  518,  518,  518,  518,
                                         518,  518,  519,  519,  519,  519,  519,  520,  520,  520,  520,  520,
                                         520,  520,  520,  520,  520,  520,  520,  520,  521,  521,  521,  521,
                                         521,  521,  521,  521,  521,  522,  522,  522,  522,  522,  522,  522,
                                         522,  522,  522,  523,  523,  523,  523,  523,  523,  523,  524,  524,
                                         524,  524,  524,  524,  524,  525,  525,  525,  525,  525,  525,  525,
                                         525,  526,  526,  526,  526,  526,  527,  527,  527,  527,  527,  527,
                                         527,  527,  527,  527,  528,  528,  528,  528,  528,  528,  528,  528,
                                         528,  528,  528,  529,  529,  529,  530,  530,  530,  530,  530,  530,
                                         530,  530,  530,  530,  531,  531,  531,  531,  531,  531,  531,  531,
                                         532,  532,  532,  532,  532,  532,  533,  533,  533,  533,  534,  534,
                                         534,  534,  535,  535,  535,  535,  535,  535,  535,  535,  535,  536,
                                         536,  536,  536,  536,  536,  536,  536,  536,  536,  537,  537,  537,
                                         537,  537,  537,  537,  537,  538,  538,  538,  538,  538,  538,  538,
                                         538,  538,  538,  538,  538,  539,  539,  539,  539,  539,  539,  539,
                                         539,  539,  539,  539,  540,  540,  540,  540,  540,  540,  540,  540,
                                         540,  540,  541,  541,  541,  541,  541,  541,  541,  541,  542,  542,
                                         542,  542,  542,  543,  543,  543,  543,  543,  543,  544,  544,  544,
                                         544,  544,  544,  544,  544,  544,  545,  545,  545,  545,  545,  545,
                                         546,  546,  546,  546,  546,  546,  546,  547,  547,  547,  547,  548,
                                         548,  548,  548,  548,  548,  548,  548,  549,  549,  549,  549,  549,
                                         549,  549,  549,  549,  549,  550,  550,  550,  550,  550,  551,  551,
                                         552,  552,  552,  552,  552,  552,  552,  552,  552,  552,  553,  553,
                                         553,  553,  553,  553,  553,  553,  554,  554,  554,  554,  554,  554,
                                         554,  555,  555,  555,  555,  555,  555,  555,  555,  555,  556,  556,
                                         556,  556,  556,  556,  556,  557,  557,  557,  557,  557,  557,  557,
                                         557,  558,  558,  558,  558,  558,  558,  558,  558,  559,  559,  559,
                                         559,  559,  559,  560,  560,  560,  560,  560,  560,  561,  561,  561,
                                         561,  561,  561,  561,  561,  561,  562,  562,  562,  562,  562,  562,
                                         562,  562,  562,  563,  563,  563,  563,  563,  563,  563,  564,  564,
                                         564,  564,  564,  564,  564,  564,  565,  565,  565,  565,  565,  566,
                                         566,  566,  566,  566,  566,  567,  567,  567,  567,  567,  567,  567,
                                         568,  568,  568,  568,  568,  569,  569,  569,  569,  569,  569,  569,
                                         570,  570,  570,  570,  570,  570,  571,  571,  571,  571,  571,  571,
                                         571,  572,  572,  572,  572,  572,  572,  572,  572,  572,  573,  573,
                                         573,  573,  573,  573,  573,  573,  573,  574,  574,  574,  574,  574,
                                         575,  575,  575,  575,  575,  575,  576,  576,  576,  576,  576,  576,
                                         576,  577,  577,  577,  577,  577,  578,  578,  578,  578,  578,  578,
                                         579,  579,  579,  579,  579,  579,  579,  579,  579,  579,  579,  579,
                                         579,  580,  580,  580,  580,  580,  580,  580,  580,  581,  581,  581,
                                         581,  581,  581,  581,  582,  582,  582,  582,  582,  582,  582,  582,
                                         582,  583,  583,  583,  583,  584,  584,  584,  584,  584,  585,  585,
                                         585,  585,  585,  586,  586,  586,  586,  586,  586,  586,  586,  587,
                                         587,  587,  587,  588,  588,  588,  588,  588,  588,  589,  589,  589,
                                         589,  589,  589,  589,  590,  590,  590,  590,  590,  590,  590,  591,
                                         591,  591,  591,  591,  591,  591,  591,  591,  592,  592,  592,  592,
                                         593,  593,  593,  593,  593,  593,  593,  593,  593,  593,  593,  594,
                                         594,  594,  594,  594,  595,  595,  595,  595,  596,  596,  596,  596,
                                         596,  596,  596,  596,  596,  596,  596,  597,  597,  597,  597,  597,
                                         597,  597,  597,  597,  597,  598,  598,  598,  598,  598,  599,  599,
                                         599,  599,  599,  599,  599,  599,  600,  600,  600,  600,  600,  601,
                                         601,  601,  601,  601,  601,  601,  601,  601,  602,  602,  602,  602,
                                         602,  602,  602,  603,  603,  603,  603,  603,  603,  604,  604,  604,
                                         604,  604,  604,  604,  605,  605,  605,  605,  605,  605,  606,  606,
                                         606,  606,  607,  607,  607,  607,  607,  607,  607,  608,  608,  608,
                                         608,  608,  608,  608,  609,  609,  609,  609,  609,  609,  609,  610,
                                         610,  610,  610,  610,  610,  611,  611,  611,  611,  611,  611,  612,
                                         612,  612,  612,  612,  612,  612,  612,  612,  612,  612,  612,  613,
                                         613,  613,  613,  613,  613,  613,  613,  614,  614,  614,  614,  614,
                                         614,  615,  615,  615,  615,  615,  616,  616,  616,  616,  617,  617,
                                         617,  617,  618,  618,  618,  618,  618,  618,  618,  618,  619,  619,
                                         619,  619,  619,  619,  619,  619,  620,  620,  620,  620,  620,  621,
                                         621,  621,  621,  621,  622,  622,  622,  622,  622,  622,  623,  623,
                                         623,  623,  623,  623,  623,  623,  623,  624,  624,  624,  624,  624,
                                         624,  624,  624,  624,  624,  625,  625,  625,  625,  625,  625,  626,
                                         626,  626,  627,  627,  627,  627,  627,  627,  628,  628,  628,  628,
                                         629,  629,  629,  629,  629,  629,  629,  630,  630,  630,  630,  630,
                                         630,  630,  630,  630,  630,  630,  631,  631,  631,  631,  631,  631,
                                         631,  632,  632,  632,  633,  633,  633,  633,  633,  633,  633,  633,
                                         634,  634,  634,  634,  635,  635,  635,  635,  636,  636,  636,  636,
                                         636,  636,  636,  637,  637,  637,  637,  637,  637,  637,  638,  638,
                                         638,  638,  638,  638,  638,  638,  638,  638,  638,  638,  639,  639,
                                         639,  639,  639,  640,  640,  640,  640,  640,  640,  641,  641,  641,
                                         641,  641,  641,  641,  641,  642,  642,  642,  642,  643,  643,  643,
                                         643,  643,  643,  643,  644,  644,  644,  644,  644,  644,  644,  644,
                                         644,  645,  645,  645,  645,  645,  645,  645,  646,  646,  646,  646,
                                         646,  646,  646,  647,  647,  647,  647,  647,  647,  648,  648,  648,
                                         648,  648,  648,  648,  649,  649,  649,  649,  649,  649,  649,  650,
                                         650,  650,  650,  650,  650,  650,  651,  651,  651,  651,  652,  652,
                                         652,  652,  652,  653,  653,  653,  653,  653,  653,  654,  654,  654,
                                         654,  654,  654,  654,  654,  655,  655,  655,  655,  655,  655,  655,
                                         655,  656,  656,  656,  656,  656,  657,  657,  657,  657,  657,  657,
                                         657,  657,  658,  658,  658,  658,  658,  658,  658,  659,  659,  659,
                                         659,  659,  660,  660,  660,  660,  660,  660,  660,  660,  660,  661,
                                         661,  661,  661,  661,  661,  661,  661,  662,  662,  662,  662,  662,
                                         662,  662,  663,  663,  663,  663,  663,  663,  663,  663,  664,  664,
                                         664,  664,  664,  664,  665,  665,  665,  665,  665,  665,  666,  666,
                                         666,  666,  666,  666,  667,  667,  667,  667,  667,  667,  667,  667,
                                         667,  668,  668,  669,  669,  669,  669,  669,  670,  670,  670,  670,
                                         670,  670,  670,  671,  671,  672,  672,  672,  672,  672,  672,  672,
                                         673,  673,  673,  673,  673,  673,  673,  673,  674,  674,  674,  674,
                                         674,  674,  674,  674,  675,  675,  675,  675,  675,  675,  675,  675,
                                         676,  676,  676,  676,  676,  676,  677,  677,  677,  677,  677,  677,
                                         678,  678,  678,  678,  678,  678,  678,  678,  678,  679,  679,  679,
                                         679,  679,  679,  680,  680,  680,  680,  681,  681,  681,  681,  681,
                                         682,  682,  682,  682,  682,  682,  682,  683,  683,  683,  683,  683,
                                         684,  684,  684,  684,  684,  684,  685,  685,  685,  685,  685,  685,
                                         686,  686,  686,  687,  687,  687,  687,  687,  687,  687,  688,  688,
                                         688,  688,  688,  688,  689,  689,  689,  689,  689,  689,  689,  689,
                                         690,  690,  690,  690,  690,  691,  691,  691,  691,  691,  692,  692,
                                         692,  692,  692,  692,  692,  692,  693,  693,  693,  693,  694,  694,
                                         694,  694,  694,  694,  694,  694,  694,  694,  694,  694,  694,  694,
                                         695,  695,  695,  696,  696,  696,  696,  696,  696,  697,  697,  697,
                                         697,  697,  697,  697,  697,  698,  698,  698,  699,  699,  699,  699,
                                         699,  699,  699,  699,  700,  700,  701,  701,  701,  701,  701,  702,
                                         702,  702,  702,  702,  702,  702,  702,  702,  702,  702,  703,  703,
                                         703,  703,  704,  704,  704,  704,  704,  704,  704,  704,  705,  705,
                                         705,  705,  705,  705,  705,  705,  706,  706,  706,  706,  706,  706,
                                         707,  707,  707,  707,  707,  707,  707,  707,  708,  708,  708,  708,
                                         708,  709,  709,  709,  709,  709,  709,  710,  710,  710,  710,  710,
                                         710,  710,  710,  711,  711,  711,  711,  711,  712,  712,  712,  712,
                                         712,  712,  712,  713,  713,  713,  713,  713,  714,  714,  714,  714,
                                         714,  715,  715,  715,  715,  716,  716,  716,  716,  717,  717,  717,
                                         717,  718,  718,  718,  718,  719,  719,  719,  719,  719,  719,  719,
                                         719,  719,  719,  720,  720,  720,  720,  721,  721,  721,  721,  721,
                                         721,  721,  721,  721,  722,  722,  722,  722,  722,  722,  722,  722,
                                         722,  723,  723,  723,  724,  724,  724,  724,  724,  724,  725,  725,
                                         725,  725,  726,  726,  726,  726,  726,  726,  726,  726,  726,  727,
                                         727,  727,  727,  727,  728,  728,  729,  729,  729,  730,  730,  730,
                                         730,  730,  730,  731,  731,  731,  732,  732,  732,  732,  732,  733,
                                         733,  733,  733,  733,  734,  734,  734,  734,  734,  735,  735,  735,
                                         735,  735,  736,  736,  736,  736,  736,  736,  736,  736,  736,  736,
                                         737,  737,  738,  738,  738,  738,  738,  738,  738,  739,  739,  739,
                                         739,  739,  739,  739,  739,  739,  739,  740,  740,  740,  741,  741,
                                         741,  741,  741,  741,  741,  741,  741,  741,  742,  742,  743,  743,
                                         743,  743,  743,  743,  743,  743,  744,  744,  744,  744,  744,  744,
                                         744,  745,  746,  746,  747,  747,  747,  748,  748,  748,  748,  748,
                                         748,  748,  749,  749,  749,  749,  750,  750,  750,  751,  751,  751,
                                         751,  751,  751,  751,  752,  752,  752,  752,  753,  753,  753,  753,
                                         753,  753,  753,  754,  754,  754,  754,  755,  755,  755,  755,  756,
                                         756,  756,  756,  756,  756,  756,  757,  757,  757,  757,  757,  757,
                                         758,  758,  758,  758,  758,  758,  759,  759,  759,  760,  760,  760,
                                         760,  760,  760,  760,  761,  761,  761,  761,  761,  761,  761,  762,
                                         762,  762,  762,  763,  763,  763,  763,  763,  764,  764,  764,  764,
                                         764,  764,  765,  765,  765,  765,  766,  766,  766,  766,  766,  766,
                                         766,  767,  767,  767,  767,  768,  768,  768,  768,  768,  768,  769,
                                         769,  769,  769,  770,  770,  770,  770,  770,  770,  770,  771,  771,
                                         771,  771,  772,  773,  773,  773,  773,  773,  773,  773,  773,  773,
                                         773,  773,  773,  773,  774,  774,  774,  774,  775,  775,  775,  775,
                                         776,  776,  776,  776,  776,  776,  776,  776,  776,  777,  777,  777,
                                         777,  778,  778,  778,  778,  779,  779,  779,  779,  779,  779,  780,
                                         780,  780,  780,  780,  781,  781,  781,  781,  781,  781,  781,  781,
                                         781,  782,  782,  782,  782,  783,  783,  784,  784,  784,  784,  784,
                                         785,  785,  785,  785,  786,  786,  786,  786,  787,  787,  787,  787,
                                         787,  787,  787,  788,  789,  789,  789,  789,  789,  790,  790,  790,
                                         790,  790,  790,  791,  791,  791,  791,  792,  792,  792,  792,  792,
                                         792,  792,  793,  793,  794,  794,  794,  794,  794,  794,  794,  795,
                                         795,  795,  795,  795,  795,  796,  796,  796,  796,  796,  796,  796,
                                         796,  796,  797,  797,  797,  797,  797,  798,  798,  798,  798,  799,
                                         799,  799,  799,  799,  800,  800,  800,  800,  800,  801,  801,  802,
                                         802,  802,  802,  802,  802,  802,  803,  803,  804,  804,  804,  804,
                                         805,  805,  805,  805,  806,  806,  806,  807,  807,  807,  807,  807,
                                         807,  807,  807,  807,  808,  808,  808,  809,  809,  809,  809,  809,
                                         809,  809,  810,  810,  810,  810,  810,  811,  811,  811,  811,  811,
                                         812,  812,  812,  813,  813,  813,  813,  813,  813,  813,  814,  814,
                                         814,  815,  815,  815,  815,  815,  816,  816,  816,  816,  816,  816,
                                         816,  816,  816,  817,  817,  817,  818,  818,  818,  818,  819,  819,
                                         819,  819,  819,  819,  819,  820,  821,  821,  821,  822,  822,  822,
                                         822,  822,  822,  823,  823,  823,  824,  824,  824,  824,  824,  824,
                                         825,  825,  825,  825,  825,  826,  826,  827,  827,  827,  827,  827,
                                         827,  827,  827,  828,  828,  828,  828,  828,  829,  829,  829,  829,
                                         830,  830,  830,  830,  830,  830,  830,  830,  831,  831,  831,  831,
                                         831,  831,  832,  832,  832,  832,  833,  833,  833,  833,  833,  834,
                                         834,  835,  835,  835,  835,  835,  835,  836,  836,  836,  836,  836,
                                         836,  837,  837,  837,  837,  837,  837,  837,  839,  839,  839,  839,
                                         839,  839,  840,  840,  840,  840,  840,  841,  841,  842,  842,  842,
                                         842,  842,  843,  843,  843,  843,  843,  844,  844,  844,  844,  844,
                                         845,  845,  845,  846,  846,  846,  846,  846,  846,  846,  846,  846,
                                         846,  847,  848,  848,  848,  848,  849,  849,  849,  850,  850,  850,
                                         850,  850,  850,  851,  851,  851,  851,  851,  852,  853,  853,  853,
                                         854,  854,  854,  854,  854,  855,  855,  855,  855,  855,  855,  855,
                                         856,  856,  857,  857,  857,  857,  858,  858,  858,  858,  858,  859,
                                         859,  859,  859,  859,  860,  860,  860,  861,  861,  861,  861,  861,
                                         862,  862,  862,  862,  863,  863,  863,  864,  864,  864,  864,  864,
                                         864,  864,  864,  865,  865,  866,  866,  866,  866,  867,  867,  867,
                                         867,  868,  868,  868,  868,  869,  869,  869,  869,  870,  870,  870,
                                         870,  871,  871,  871,  871,  871,  871,  871,  872,  872,  873,  873,
                                         873,  873,  873,  874,  874,  874,  875,  875,  875,  875,  876,  876,
                                         876,  876,  876,  876,  876,  876,  877,  877,  877,  877,  877,  877,
                                         877,  878,  878,  878,  878,  878,  878,  878,  879,  879,  879,  879,
                                         880,  880,  880,  880,  880,  881,  881,  881,  882,  882,  882,  882,
                                         883,  883,  884,  884,  885,  885,  885,  885,  885,  886,  886,  886,
                                         887,  887,  887,  888,  888,  888,  888,  889,  889,  890,  890,  890,
                                         891,  891,  891,  891,  891,  891,  891,  891,  892,  892,  893,  893,
                                         893,  893,  893,  894,  894,  894,  894,  894,  895,  895,  895,  895,
                                         895,  896,  896,  896,  896,  896,  896,  896,  896,  897,  897,  898,
                                         898,  898,  898,  899,  899,  899,  900,  900,  901,  901,  901,  903,
                                         904,  904,  904,  905,  905,  905,  905,  906,  906,  906,  906,  907,
                                         907,  907,  907,  908,  908,  908,  909,  909,  909,  909,  909,  909,
                                         909,  910,  910,  910,  910,  910,  910,  910,  911,  911,  911,  912,
                                         912,  912,  913,  913,  914,  914,  914,  914,  914,  914,  914,  915,
                                         915,  916,  916,  916,  916,  916,  916,  917,  917,  918,  918,  920,
                                         920,  920,  921,  921,  921,  921,  921,  923,  923,  923,  923,  923,
                                         924,  924,  924,  925,  925,  925,  925,  925,  925,  926,  926,  926,
                                         926,  926,  926,  927,  927,  927,  928,  928,  928,  929,  929,  929,
                                         929,  929,  930,  930,  930,  930,  931,  931,  932,  932,  932,  933,
                                         933,  933,  933,  934,  934,  934,  935,  935,  935,  936,  936,  936,
                                         937,  937,  937,  937,  938,  939,  939,  939,  940,  940,  940,  940,
                                         940,  940,  941,  941,  942,  942,  942,  942,  942,  943,  943,  944,
                                         944,  944,  945,  945,  945,  945,  946,  946,  946,  946,  947,  947,
                                         947,  948,  948,  948,  948,  949,  949,  950,  950,  950,  951,  951,
                                         951,  951,  951,  951,  952,  952,  952,  953,  953,  955,  955,  955,
                                         955,  955,  956,  956,  957,  957,  957,  957,  958,  958,  958,  959,
                                         959,  960,  960,  960,  960,  961,  961,  962,  962,  962,  963,  963,
                                         963,  963,  964,  964,  964,  964,  964,  965,  965,  965,  965,  965,
                                         965,  965,  966,  966,  966,  967,  967,  968,  968,  968,  969,  969,
                                         969,  969,  970,  970,  970,  970,  971,  971,  971,  972,  972,  973,
                                         973,  973,  973,  974,  974,  974,  975,  975,  975,  976,  976,  976,
                                         977,  977,  978,  978,  978,  979,  979,  979,  979,  980,  980,  980,
                                         980,  980,  981,  981,  981,  981,  982,  982,  983,  983,  983,  984,
                                         984,  985,  985,  985,  985,  985,  985,  986,  986,  987,  988,  988,
                                         988,  988,  989,  989,  990,  990,  991,  991,  991,  992,  992,  992,
                                         992,  993,  993,  993,  993,  994,  995,  995,  996,  996,  996,  996,
                                         996,  996,  997,  997,  997,  997,  998,  998,  999,  999,  999,  999,
                                         999,  999,  999,  999, 1000, 1001, 1001, 1002, 1002, 1002, 1002, 1002,
                                     1003, 1003, 1004, 1004, 1004, 1005, 1005, 1005, 1006, 1006, 1006, 1007,
                                     1007, 1007, 1007, 1008, 1008, 1009, 1009, 1010, 1010, 1010, 1010, 1010,
                                     1011, 1011, 1011, 1012, 1012, 1012, 1012, 1014, 1014, 1014, 1014, 1015,
                                     1015, 1015, 1015, 1016, 1016, 1016, 1017, 1017, 1017, 1017, 1017, 1018,
                                     1018, 1019, 1019, 1019, 1020, 1020, 1021, 1021, 1021, 1021, 1021, 1022,
                                     1023, 1023, 1024, 1024, 1024, 1024, 1024, 1025, 1026, 1026, 1027, 1027,
                                     1027, 1027, 1028, 1028, 1028, 1028, 1029, 1029, 1029, 1029, 1029, 1030,
                                     1030, 1030, 1030, 1031, 1031, 1032, 1032, 1033, 1034, 1034, 1034, 1034,
                                     1035, 1035, 1035, 1035, 1035, 1036, 1036, 1036, 1036, 1038, 1039, 1039,
                                     1039, 1039, 1039, 1040, 1041, 1041, 1042, 1042, 1042, 1043, 1043, 1043,
                                     1043, 1043, 1043, 1044, 1044, 1045, 1045, 1046, 1046, 1046, 1046, 1047,
                                     1047, 1048, 1048, 1049, 1049, 1050, 1050, 1051, 1051, 1051, 1051, 1051,
                                     1052, 1052, 1052, 1053, 1053, 1054, 1054, 1054, 1054, 1055, 1055, 1056,
                                     1057, 1057, 1057, 1058, 1058, 1058, 1058, 1058, 1058, 1060, 1060, 1060,
                                     1060, 1061, 1061, 1061, 1062, 1062, 1062, 1063, 1063, 1065, 1065, 1065,
                                     1065, 1066, 1067, 1067, 1067, 1068, 1068, 1069, 1069, 1069, 1070, 1070,
                                     1070, 1070, 1070, 1070, 1072, 1073, 1073, 1073, 1073, 1074, 1075, 1075,
                                     1075, 1076, 1076, 1076, 1076, 1076, 1077, 1077, 1078, 1078, 1079, 1079,
                                     1079, 1080, 1080, 1081, 1081, 1082, 1084, 1084, 1084, 1085, 1085, 1086,
                                     1087, 1087, 1087, 1088, 1088, 1088, 1089, 1090, 1090, 1090, 1091, 1091,
                                     1092, 1093, 1093, 1093, 1094, 1094, 1094, 1094, 1094, 1095, 1095, 1096,
                                     1097, 1097, 1098, 1098, 1098, 1099, 1099, 1100, 1100, 1100, 1101, 1101,
                                     1102, 1102, 1103, 1104, 1104, 1105, 1105, 1105, 1106, 1107, 1107, 1107,
                                     1108, 1108, 1108, 1109, 1109, 1110, 1110, 1110, 1111, 1111, 1112, 1113,
                                     1113, 1113, 1115, 1116, 1116, 1116, 1117, 1117, 1117, 1119, 1119, 1121,
                                     1122, 1124, 1124, 1124, 1124, 1124, 1124, 1125, 1127, 1127, 1127, 1128,
                                     1128, 1129, 1130, 1130, 1130, 1131, 1131, 1131, 1132, 1132, 1133, 1134,
                                     1134, 1134, 1134, 1135, 1136, 1136, 1137, 1137, 1138, 1139, 1139, 1139,
                                     1139, 1140, 1140, 1141, 1142, 1142, 1142, 1143, 1143, 1144, 1144, 1145,
                                     1145, 1146, 1146, 1147, 1147, 1148, 1149, 1149, 1149, 1149, 1150, 1150,
                                     1151, 1151, 1152, 1153, 1154, 1155, 1156, 1156, 1157, 1158, 1158, 1158,
                                     1158, 1161, 1161, 1163, 1163, 1163, 1164, 1164, 1164, 1165, 1165, 1165,
                                     1165, 1167, 1167, 1168, 1168, 1168, 1170, 1171, 1171, 1171, 1172, 1173,
                                     1174, 1174, 1174, 1176, 1176, 1176, 1176, 1178, 1179, 1180, 1181, 1182,
                                     1182, 1183, 1183, 1183, 1185, 1185, 1185, 1185, 1187, 1188, 1188, 1189,
                                     1190, 1190, 1191, 1192, 1192, 1192, 1193, 1195, 1195, 1196, 1198, 1198,
                                     1198, 1199, 1199, 1201, 1201, 1201, 1202, 1202, 1203, 1204, 1205, 1207,
                                     1209, 1210, 1212, 1212, 1215, 1216, 1218, 1219, 1219, 1219, 1222, 1222,
                                     1222, 1222, 1225, 1227, 1228, 1229, 1230, 1231, 1231, 1233, 1234, 1236,
                                     1238, 1239, 1245, 1246, 1249, 1253, 1256, 1258, 1265
                                       ]))

        def test_ba11a(self):
            '''BA11A Construct the Graph of a Spectrum'''
            graph = SpectrumGraph([57, 71, 154, 185, 301, 332, 415, 429, 486])
            self.assertEqual(9,len(graph))
        @skip('')
        def test_ba11b(self):
            '''BA11B Implement DecodingIdealSpectrum'''
            self.assertEqual('GPFNA',DecodeIdealSpectrum([57, 71, 154, 185, 301, 332, 415, 429, 486]))

        @skip('Code up test')
        def test_ba11c(self):
            '''BA11C Convert a Peptide into a Peptide Vector'''


        def test_ba11d(self):
            '''BA11D Convert a Peptide Vector into a Peptide'''
            self.assertAlmostEqual('XZZXX',CreatePeptide([0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]))

        @skip('')
        def test_ba11e(self):
            '''BA11E Sequence a Peptide'''
            self.assertEqual('XZZXX',
                             SequencePeptide([0, 0, 0, 4, -2, -3, -1, -7, 6, 5, 3, 2, 1, 9, 3, -8, 0, 3, 1, 2, 1, 0],
                                             protein_masses=test_masses))

        def test_splc(self):
            '''SPLC	RNA Splicing'''
            string='''>Rosalind_10
            ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
            >Rosalind_12
            ATCGGTCGAA
            >Rosalind_15
            ATCGGTCGAGCGTGT'''
            self.assertEqual('MVYIADKQHVASREAYGHMFKVCA',splc(FastaContent(string.split('\n'))))

        def test_prot(self):
            ''' PROT Translating RNA into Protein'''
            self.assertEqual('MAMAPRTEINSTRING',prot('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'))

    main()
