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

# DecodeIdealSpectrum
#
# Reconstruct a peptide from its ideal spectrum.
#
# Input: A  list of integers, Spectrum.
#
# Return: An amino acid string with an ideal spectrum that matches Spectrum.

def DecodeIdealSpectrum(Spectrum):
    def bfs(adj,paths = [[0]]):
        while True:
            new_paths=[]
            for path in paths:
                a = path[-1]
                if a in adj:
                    for b,_ in adj[a]:
                        new_paths.append(path + [b])
            if len(new_paths)==0:
                return paths
            else:
                paths = new_paths

    adj = SpectrumGraph(Spectrum)
    paths = bfs(adj)
    return ''