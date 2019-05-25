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
    masses = [extended_masses[p] for p in peptide]
    result = []
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

def conv(S,T,eps=0.001):
    diff = sorted([s-t for s in S for t in T])
    accumulated = []
    latest      = diff[0]
    count       = 1
    
    for term in diff[1:]+[666]:
        if abs(term-latest)<eps:
            count+=1
        else:
            accumulated.append((count,latest))
            latest = term
            count  = 1
    index_max = argmax([i for i,_ in accumulated])
    c,v = accumulated[index_max]
    return c,v

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