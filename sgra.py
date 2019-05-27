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
        #3524.8542,
        #3623.5245,
        #3710.9335,
        #3841.974,
        #3929.00603,
        #3970.0326,
        #4026.05879,
        #4057.0646,
        #4083.08025
        3895.22521511,
        3966.26232511,
        4117.87918874,
        4134.04602924,
        4152.34163511,
        4180.09378653,
        4265.42569511,
        4281.99952562,
        4326.09226024,
        4326.38297797,
        4428.48902511,
        4442.90766721,
        4614.56833511,
        4618.2773078,
        4622.99505648,
        4704.88889633,
        4727.65239511,
        4798.02914636,
        4824.70515511,
        4952.76373511,
        5049.81649511,
        5052.07485394,
        5058.73230054,
        5106.4317037,
        5106.83795511,
        5178.0444579,
        5234.0817192,
        5243.89686511,
        5314.93397511,
        5414.00238511,
        5456.12488616,
        5469.45213693,
        5551.06129511,
        5626.20477754,
        5666.08823511,
        5689.1975723,
        5735.94438895,
        5744.3584684,
        5781.54132697,
        5794.14681511,
        5803.93213381,
        5850.14916782,
        5851.16827511,
        5908.8298212,
        5980.21086511,
        5980.21086511,
        6093.29492511,
        6117.26977511,
        6135.62286963,
        6222.8259566,
        6256.35825511,
        6303.34908511,
        6326.44069466,
        6387.39874511,
        6418.37602511,
        6501.44167511,
        6562.57403961,
        6574.47713511,
        6594.64048783,
        6630.48426511,
        6631.49859511,
        6640.50861501,
        6698.16337523,
        6731.53194511,
        6734.50778511,
        6835.55546511,
        6917.61125511,
        6950.58240511,
        7014.66401511,
        7049.65081511,
        7115.71169511,
        7162.73487511,
        7243.77027511,
        7332.03722191,
        7374.81076511,
        7391.86847479,
        7459.81271703,
        7505.85125511,
        7562.87271511,
        7676.91564511,
        7699.9610851,
        7805.01060511,
        7818.3419703,
        7961.11171511,
        8058.16447511,
        8205.23288511,
        8290.83449049,
        8334.27547511
        
                 ]))