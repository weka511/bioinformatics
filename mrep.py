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

#  MREP Identifying Maximal Repeats

import argparse
import os
import time
from   helpers import read_strings

class SuffixNode:
    def __init__(self,suffix_link=None):
        self.children = {}
        if suffix_link==None:
            self.suffix_link = self
        else:
            self.suffix_link = suffix_link
    def add_link(self,c,v):
        self.children[c] = v    
        
    def search(self,
               path        = [],
               visit       = lambda node,path: print (f'{path}:{node}'),
               extend_path = lambda node,path,child: path + [child]):
        print ( '**',self.children)
        for child in self.children:
            visit(child,path)
            path1 = extend_path(self,path,child)
            self.children[child].search(path1,extend_path=extend_path)
            
# build_suffix_trie
        
def build_suffix_trie(s):
    assert len(s)>0
    Root    = SuffixNode()
    Longest = SuffixNode(suffix_link=Root)
    Root.add_link(s[0],Longest)
    for c in s[1:]:
        Current  = Longest
        Previous = None
        while c not in Current.children:
            r1 = SuffixNode()
            Current.add_link(c,r1)
            if Previous is not None:
                Previous.suffix_link = r1
            Previous = r1
            Current  = Current.suffix_link
        Previous.suffix_link = Root if Current is Root else  Current.children[c]
        Longest = Longest.children[c]
    return Root
    
    
def mrep(s,minimum_length=20):
    def create_candidates(k):
        Candidates = []
        for i in range(len(Suffixes)):
            s1,_    = Suffixes[i]
            if len(s1)<k: continue
            if len(Candidates)==0:
                Candidates.append([i])
            else:
                indices = Candidates[-1]
                s2,_    = Suffixes[indices[0]]
                if s1[0:k]== s2[0:k]:
                    Candidates[-1].append(i)
                else:
                    Candidates.append([i])
        return [indices for indices in Candidates if len(indices)>1]
        #return [(k,indices) for indices in Candidates if len(indices)>1]
    def matches(k,index,target):
        suffix,_ = Suffixes[index] 
        return k<len(suffix) and k<len(target) and suffix[k]==target[k]
    
    Suffixes   = sorted([(s[i:],i) for i in range(len(s))])
    Candidates = create_candidates(minimum_length)

    for k in range(minimum_length,len(s)-minimum_length):
        Extended = []
        for Indices in Candidates:
            target,_ = Suffixes[Indices[0]]
            if all(matches(k,index,target) for index in Indices[1:]):
                Extended.append(Indices)
            else:
                s1,_ = Suffixes[Indices[0]]
                yield s1[0:k]
        Candidates = [indices for indices in Extended]


if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('MREP Identifying Maximal Repeats')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    
    if args.sample:
        
        #for r in mrep('TAGTTAGCGAGA',minimum_length=2):
            #print (r)
        for r in mrep('TAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTATTATATAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTAT'):
            print (r)
        
    

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = mrep(Input[0])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
