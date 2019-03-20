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

from reference_tables import integer_masses
from graphs import dfs

# SpectrumGraph
#
# Construct the Graph of a Spectrum
# Input: A space-delimited list of integers Spectrum.
#
# Return: Graph(Spectrum)

def SpectrumGraph(spectrum):
    # add
    #
    # Add one point to graph
    
    def add(index=-1):
        value = spectrum[index] if index>-1 else 0
        for j in range(index+1,len(spectrum)):
            diff = spectrum[j]-value
            if diff in reversed:
                if not value in product:
                    product[value]=[]
                for protein in reversed[diff]:
                    product[value].append((spectrum[j],protein))
    reversed={}
    for k,v in integer_masses.items():
        if not v in reversed:
            reversed[v]=[k]
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

def CreatePeptideVector(peptide):
    extended_masses={'X':4,'Z':5}
    extended_masses.update(integer_masses)
    masses = [extended_masses[p] for p in peptide]
    result = []
    for m in masses:
        result = result + ([0]*(m-1))
        result.append(1)
    return result