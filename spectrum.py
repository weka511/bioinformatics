#    Copyright (C) 2019 Greenweaves Software Limited
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
#
#    Utilities for mass sprctroscopy

from reference_tables import integer_masses,amino_acids
from graphs import dfs
from bisect import bisect
from helpers import create_list
from numpy import argmax

# SpectrumGraph
#
# Construct the Graph of a Spectrum
# Input: A space-delimited list of integers Spectrum.
#
# Return: Graph(Spectrum)

def invert(masses):
    inverted={}
    for k,v in masses.items():
        if not v in inverted:
            inverted[v]=[k]
    return inverted

def SpectrumGraph(spectrum):
    # add
    #
    # Add one point to graph
    inverted = invert(integer_masses)
    def add(index=-1):
        value = spectrum[index] if index>-1 else 0
        for j in range(index+1,len(spectrum)):
            diff = spectrum[j]-value
            if diff in inverted:
                if not value in product:
                    product[value]=[]
                for protein in inverted[diff]:
                    product[value].append((spectrum[j],protein))

    product = {}
    add()
    for i in range(len(spectrum)):
        add(i)
    return product

def linearSpectrum(peptide):
    def getSpectrum(peptide):
        spectrum = set()
        for i in range(len(peptide)):
            spectrum.add(sum(peptide[:i]))
            spectrum.add(sum(peptide[i:]))
     
        return sorted(list(spectrum))
    return getSpectrum( [integer_masses[p] for p in peptide]) if peptide.isalpha() else getSpectrum(peptide)

# DecodeIdealSpectrum
#
# Reconstruct a peptide from its ideal spectrum.
#
# Input: A  list of integers, Spectrum.
#
# Return: An amino acid string with an ideal spectrum that matches Spectrum.

def DecodeIdealSpectrum(Spectrum):
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
        peptide            = ''.join(s for (_,s) in path)
        generated_spectrum = linearSpectrum(peptide)
        if generated_spectrum[1:]==Spectrum:
            return peptide
 
    return None

# CreatePeptideVector
#
# Convert a Peptide into a Peptide Vector 
#
# Convert a peptide into a binary peptide vector.
#
# Given an amino acid string Peptide = a1 . . . an of length n, we will represent its
# prefix masses using a binary peptide vector Peptide' with mass(Peptide) coordinates.
# This vector contains a 1 at each of the n prefix coordinates
#
# mass(a1), mass(a1 a2), . . . , mass(a1 a2 . . . an ) ,
# and it contains a 0 in each of the remaining noise coordinates.
#
# Input: A peptide P.
#
# Return: The peptide vector of P.
#
# Note: In this chapter, all dataset problems implicitly use the standard integer-valued mass
# table for the regular twenty amino acids. Examples sometimes use imaginary amino
# acids X and Z having respective integer masses 4 and 5.

def create_extended():
    extended_masses = {'X':4,'Z':5}
    extended_masses.update(integer_masses)
    return extended_masses

def CreatePeptideVector(peptide):
    extended_masses = create_extended()
    masses          = [extended_masses[p] for p in peptide]
    result          = []
    for m in masses:
        result = result + ([0]*(m-1))
        result.append(1)
    return result

# CreatePeptide
#
# Convert a Peptide Vector into a Peptide
#
def CreatePeptide(vector):
    extended_masses = create_extended()
    masses_offset   = [i+1 for i in range(len(vector)) if vector[i]>0]
    masses          = [b-a for (a,b) in zip([0]+masses_offset[:-1],masses_offset)]
    inverted_masses = invert(extended_masses)
    return ''.join( [str(inverted_masses[m][0]) for m in masses])

# conv
#
# Comparing Spectra with the Spectral Convolution
#
# Comparing Spectraclick to collapse
#
# Suppose you have two mass spectra, and you want to check if they both were obtained from the same protein; 
# you will need some notion of spectra similarity. The simplest possible metric would be to count the number
# of peaks in the mass spectrum that the spectra share, called the shared peaks count; 
# its analogue for simplified spectra is the number of masses that the two spectra have in common.

# The shared peaks count can be useful in the simplest cases, but it does not help us if, for example,
# one spectrum corresponds to a peptide contained inside of another peptide from which the second
# spectrum was obtained. In this case, the two spectra are very similar, but the shared peaks count
# will be very small. However, if we shift one spectrum to the right or left, then shared peaks will align.
# In the case of simplified spectra, this means that there is some shift value `x` such that adding
# x to the weight of every element in one spectrum should create a large number of matches in the other spectrum.
#
# Inputs: Two multisets of positive real numbers S1 and S2
#
# The size of each multiset is at most 200.

# Return: The largest multiplicity of S1- S2, as well as the absolute value of the number x 
# maximizing (S1-S2)(x) (you may return any such value if multiple solutions exist).

def conv(S,T,eps=0.001):
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
            
    return accumulated[argmax([i for i,_ in accumulated])]
    

def create_lookup(amino_acids=amino_acids):
    pairs = sorted([(abbrev,value.mon_mass) for abbrev,value in amino_acids.items()],
                   key =lambda x:x[1])
    pairs.append(('?',999999999999))
    masses = [mass for (_,mass) in pairs]
    return masses,pairs    

def get_abbrev(diff,masses,pairs):
    index = bisect(masses,diff)
    m1 = masses[index]
    m0 = masses[(index-1) if index>0 else 0]
    if index>0 and diff-m0 < m1-diff:
        index-=1
    abbrev,_ = pairs[index]

    return abbrev

# spectrum2protein
#
# spec Inferring Protein from Spectrum
#
#    Introduction to Mass Spectrometry
#
#    In “Calculating Protein Mass”, we briefly mentioned an analytic chemical method called mass spectrometry,
#    which aims to measure the mass-to-charge ratio of a particle or a molecule. In a mass spectrometer,
#    a sample is vaporized (turned into gas), and then particles from the sample are ionized. The resulting
#    ions are placed into an electromagnetic field, which separates them based on their charge and mass.
#    The output of the mass spectrometer is a mass spectrum, or a plot of ions' possible mass-to-charge ratio 
#    values with the intensity (actual observed frequency) of ions having these mass-to-charge values.
#
#    For the moment, we will ignore charge and consider a list of the ions' monoisotopic masses as a
#    simplified spectrum. Researchers do not possess cheap technology to go in and examine a protein one
#    amino acid at a time (molecules are too submicroscopic). Instead, to determine a protein's structure,
#    we will split several copies of the protein into smaller pieces, then weigh the resulting fragments. 
#    To do this, we assume that each cut (breakage point) occurs between two amino acids and that we can 
#    measure the mass of the resulting pieces for all possible cuts.
#
#    For example, the (unknown) protein "PRTEIN" can be cut in five possible ways:
#    "P" and "RTEIN"; "PR" and "TEIN"; "PRT" and "EIN"; "PRTE" and "IN"; "PRTEI" and "N".
#    We then can measure the masses of all fragments, including the entire string. The "left" end of
#    a protein is called its N-terminus, and the ions corresponding to the protein string's prefixes 
#    (P, PR, PRT, PRTE, PRTEI) are called b-ions. The "right" end of the protein is called its C-terminus,
#    and the ions corresponding to the string's suffixes (N, IN, EIN, TEIN, RTEIN) are called y-ions. 
#    The difference in the masses of two adjacent b-ions (or y-ions) gives the mass of one amino acid residue; 
#    for example, the difference between the masses of "PRT" and "PR" must be the mass of "T." By extension,
#    knowing the masses of every b-ion of a protein allows us to deduce the protein's identity.
#
# The prefix spectrum of a weighted string is the collection of all its prefix weights.
#
# Input: A list L of n (n<=100) positive real numbers.
#
# Return: A protein string of length n-1 whose prefix spectrum is equal to L (if multiple solutions exist,
# you may output any one of them). Consult the monoisotopic mass table.

def spectrum2protein(ms):
    masses,pairs = create_lookup() 
    return ''.join([get_abbrev(diff,masses,pairs) for diff in [m1-m0 for m0,m1 in zip(ms[:-1],ms[1:])]])   

def complete_spectrum(P):
    def spectrum(S):
        return sum([amino_acids[s].mon_mass for s in S])
    prefixes = [P[:i] for i in range(1,len(P))] +[P]
    suffixes = [P[i:] for i in range(1,len(P))]
    ss= [spectrum(p) for p in prefixes + suffixes]
    return ss


#  prsm
#
# Match a Spectrum to a Protein
#
# Searching the Protein Database
#
# Many proteins have already been identified for a wide variety of organisms. Accordingly,
# there are a large number of protein databases available, and so the first step after 
# creating a mass spectrum for an unidentified protein is to search through these databases
# for a known protein with a highly similar spectrum. In this manner, many similar proteins 
# found in different species have been identified, which aids researchers in determining protein function.
#
# In “Comparing Spectra with the Spectral Convolution”, we introduced the spectral convolution 
# and used it to measure the similarity of simplified spectra. In this problem, we would like
# to extend this idea to find the most similar protein in a database to a spectrum taken from
# an unknown protein. Our plan is to use the spectral convolution to find the largest possible
# number of masses that each database protein shares with our candidate protein after shifting,
# and then select the database protein having the largest such number of shared masses.

# Inputs: A positive integer n followed by a collection of n protein strings s1, s2, ..., sn and a multiset R
# of positive numbers (corresponding to the complete spectrum of some unknown protein string).

# Return: The maximum multiplicity of R-S[sk]
# taken over all strings sk, followed by the string sk for which this 
# maximum multiplicity occurs (you may output any such value if multiple solutions exist).

def prsm(s,R):
    def count(c):
        return c[1][0]
    m = 0
    i = 0
    Ss = [(P,complete_spectrum(P)) for P in s]
 
    Cs = [(P,conv(R,S1,eps=0.00001)) for (P,S1) in Ss]

    index = argmax([count(c) for c in Cs])

    _,(b,_) = Cs[index]
    
    return b,s[index]