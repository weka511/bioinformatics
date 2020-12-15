#   Copyright (C) 2020 Greenweaves Software Limited

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

#  ALPH  Alignment-Based Phylogeny 

import argparse
import os
import time
from helpers import read_strings
from newick  import newick_to_adjacency_list
from numpy   import argmin
from fasta   import FastaContent, fasta_out

# alph
#
# Given: A rooted binary tree T on n  species, given in Newick format, followed by a multiple alignment of m 
#        augmented DNA strings having the same length (at most 300 bp) corresponding to the species
#        and given in FASTA format.
#
# Return: The minimum possible value of dH(T), followed by a collection of DNA strings to be 
#         assigned to the internal nodes of T that will minimize dH(T).

def alph(T,Alignment,Alphabet=['A','T','C','G','-']):
    
    # create_fixed_alignments
    #
    # Extract dictionary of leaves from Alignment
    #
    # Returns: length of any string in alignment, plus dictionary of leaves
    def create_fixed_alignments():
        Leaves = {}
        k      = None
        for i in range(0,len(Alignment),2):
            Leaves[Alignment[i]] = Alignment[i+1]
            if k==None:
                k = len(Alignment[i+1])
            else:
                assert k==len(Alignment[i+1]),f'Alignments should all have same length.'
                
        return k,Leaves
    
    # SmallParsimony
    #
    # This is the Small Parsimony algorihm from Pevzner and Compeau, which
    # processes a single character
    #
    # Parameters:
    #      l       Index of character in Alignment
    # Returns:     Score of best assignment, plus an assignment of character that provides this score
    
    def SmallParsimony(l):
    
        # is_ripe
        #
        # Determine whether now is ready for processing
        # A ripe node is one that hasn't been processed,
        # but its children have
        
        def is_ripe(v):
            for child in Adj[v]:
                if not Tag[child]: return False
            return True 
        
        # find_ripe
        #
        # Find list of nodes that are ready to be processed
        #
        # Input:   A list of nodes
        # Returns: Two lists, those ready for processing, and those which are not
        
        def find_ripe(Nodes):
            Ripe   = []
            Unripe = []
            for v in Nodes:
                if is_ripe(v):
                    Ripe.append(v)
                else:
                    Unripe.append(v)
            return Ripe,Unripe
        
        # get_distance
        #
        # Get total distance of node from its children assuming one trial assignmnet
        #
        # Parameters:
        #     v       Current node
        #     k       Index of character for trial
        
        def get_distance(v,k):
            # delta
            #
            # The delta function from Pevzner and Compeau: not the Kronecker delta
            
            def delta(i,j):
                return 0 if i==j else 1 
            
            # best_alignment
            #
            # Find best alignment with child (measured by varying child's index) given 
            # the current choice of character in this node
            #
            # Parameters:
            #      k       Trial alignmnet for this node
            
            def best_alignment(child):
                return min([s[child][i] + delta(i,k)  for i in range(len(Alphabet))])
            
            return sum([best_alignment(child) for child in Adj[v]])
 
        # backtrack
        #
        # Perform a depth first search through all nodes to determive alignmant. 
        # Parameters:
        #     root    Root node
        #     s       Scores for all possible best assignments to all nodes   
        # Returns: 
        #     score    Score of best assignment,
        #     ks       For each node the assignment of character that provides this score
        #              represented an an index into alphabet
        
        def backtrack(root,s):
            score = None        # Final score. Since the scores in s are cumulative
                                # we can determine this from root node
            Open  = [root]      # List of nodes that are ready to be processed
            ks    = {}          # Final assignments of characters
            while len(Open)>0:
                v = Open.pop(0)
                # Find index of best assignmnet from current node
                index = argmin([s[v][k] for k in range(len(Alphabet))])
                if score==None:
                    score = s[v][index]
                ks[v] = index
                for child in Adj[v]:
                    Open.append(child)  
                    
            return score,ks
        
        s      = {}           # Scores for nodes
        Tag    = {}           # Nodes that have been processed
        ToBeProcessed   = []           # Nodes that have yet to be processed
        
        # Partition nodes into two groups: leaves are easily processed,
        # the others are all marked as unprocessed
        for v in Adj.keys():
            if v in Leaves:
                char      = Leaves[v][l]
                s[v]      = [0 if Alphabet[k]==char else float('inf') for k in range(len(Alphabet))]
                Tag[v]    = True
            else:
                Tag[v] = False
                ToBeProcessed.append(v)
                
        Ripe,ToBeProcessed = find_ripe(ToBeProcessed)
        while len(Ripe)>0:
            for v in Ripe:
                s[v]      = [get_distance(v,k) for k in range(len(Alphabet))]
                Tag[v]    = True
            Ripe,ToBeProcessed = find_ripe(ToBeProcessed)
            
        assert len(ToBeProcessed)==0,'If there are no ripe nodes, ToBeProcessed should be exhausted'
        return backtrack(v,s)
           
    Adj      = newick_to_adjacency_list(T)   
    L,Leaves = create_fixed_alignments()
    assert len([node for node,value in Adj.items() if len(value)==0 and node not in Leaves])==0,\
           f'Some nodes have no children, but have no strings in alignment'
    
    Assignment = {a:[] for a in Adj.keys()}
    
    d = 0
    for l in range(L):
        score,ks = SmallParsimony(l)
        d       += score
        for v,index in ks.items():
            Assignment[v].append(Alphabet[index])
              
    return d,[(f'{a}',''.join(b)) for a,b in Assignment.items() if len(Adj[a])!=0]
             
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('ALPH  Alignment-Based Phylogeny')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        fc = FastaContent(['>ostrich',
                           'AC',
                           '>cat',
                           'CA',
                           '>duck',
                           'T-',
                           '>fly',
                           'GC',
                           '>elephant',
                           '-T',
                           '>pikachu',
                           'AA'
                           ])
        
        d,Assignment = alph('(((ostrich,cat)rat,(duck,fly)mouse)dog,(elephant,pikachu)hamster)robot;',fc.to_list())
        print (d)
        for label,String in Assignment:
            for line in fasta_out(label,String):
                print (line)        
        
  
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        fc     = FastaContent(Input[1:])         
        d,Assignment = alph(Input[0],fc.to_list())
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            print (d)
            f.write(f'{d}\n')
            for label,String in Assignment:
                for line in fasta_out(label,String):
                    print (line)
                    f.write(f'{line}\n')            
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
 