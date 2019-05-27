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
#    sgra Using the Spectrum Graph to Infer Peptides

from reference_tables import amino_acids
from spectrum import create_lookup,get_abbrev
from numpy import argmax

# sgra
#
# In this problem, we say that a weighted string s=s1s2...sn
# matches L if there is some increasing sequence of positive real numbers (w1,w2,,,,,wn+1)
# in L such that w(s1)=w2-w1, w(s2)=w3-w2, ..., and w(sn)=wn+1-wn
#
# Input: A list L (of length at most 100) containing positive real numbers.
#
# Return: The longest protein string that matches the spectrum graph of L (if multiple solutions exist,
#  you may output any one of them). Consult the monoisotopic mass table.
#
# NB: we want the longest path trough the spectrum graph
#
def sgra(L,Alphabet=amino_acids,epsilon=0.001):
    
    # create_spectrum_graph
    #
    # For a weighted alphabet A and a collection L of positive real numbers, the spectrum graph of L
    # is a digraph constructed in the following way. First, create a node for every real number in L.
    # Then, connect a pair of nodes with a directed edge (u,v) if v>u and v-u is equal to the weight of a single symbol in A
    #
    # We may then label the edge with this symbol.
    
    def create_spectrum_graph():
        masses,pairs   = create_lookup()
        G = {}
        for u in L:
            for v in L:
                if u<v:
                    abbrev = get_abbrev(v-u,masses,pairs)
                    if abs(v-u-amino_acids[abbrev].mon_mass)<epsilon:
                        if not u in G:
                            G[u]=[]
                        G[u].append((abbrev,v))            
        return G
    
    Outs = set()
    Runs = []
    def dfs(key,G,path):
        Runs.append(path)
        if key in Outs: return
        if not key in G: return
        for amino_acid,mass in G[key]:
            dfs(mass,G,path+[amino_acid])
            Outs.add(mass)
            
    G    = create_spectrum_graph()
    
    for key in sorted(G.keys()):
        #print (key,G[key])
        dfs(key,G,[])
 
    index           = argmax([len(run) for run in Runs])
    return ''.join(Runs[index])        

if __name__=='__main__':
    print (sgra([
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