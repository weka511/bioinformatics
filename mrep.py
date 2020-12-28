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
from   snp     import SuffixArray
  
    
def mrep(s,minimum_length=20):

    def create_candidates(k):
        Candidates = []
        for i in range(len(suffix_array)):
            j = suffix_array[i]
            if n-j<k: continue
 
            s1      = s[j:]
            if len(Candidates)==0:
                Candidates.append([i])
            else:
                indices = Candidates[-1]
                j2      = suffix_array[indices[0]]
                s2      = s[j2:]
                if s1[0:k]== s2[0:k]:
                    Candidates[-1].append(i)
                else:
                    Candidates.append([i])
        return [indices for indices in Candidates if len(indices)>1]
        
    def matches(k,index,target):
        j2       = suffix_array[index]
        suffix   = s[j2:]
        return k<len(suffix) and k<len(target) and suffix[k]==target[k]
 
    def try_right_extension(Candidates):
        Strings = []
        for k in range(minimum_length,len(s)-minimum_length):
            Extended = []
            for Indices in Candidates:
                j2      = suffix_array[Indices[0]]
                target  = s[j2:]
                if all(matches(k,index,target) for index in Indices[1:]):
                    Extended.append(Indices)
                else:
                    j2  = suffix_array[Indices[0]]
                    s1  = s[j2:]
                    Strings.append(s1[0:k])
            Candidates = [indices for indices in Extended]
        return Strings
    
    n            = len(s+'$')
    suffix_array = SuffixArray(s+'$')
    Candidates   = create_candidates(minimum_length)
    Candidates1  = try_right_extension(Candidates)
    Candidates2 = []
    for i in range(len(Candidates1)):
        subsumed = False
        for j in range(i+1,len(Candidates1)):
            if len(Candidates1[i])==len(Candidates1[j]): continue
            subsumed = Candidates1[j].find(Candidates1[i]) > -1
            if subsumed: break
        if not subsumed:
            Candidates2.append(Candidates1[i])
            
    return Candidates2


if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('MREP Identifying Maximal Repeats')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    
    if args.extra:
        for r in mrep('TAGTTAGCGAGA',minimum_length=2):
            print (r)    
            
    if args.sample:
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in mrep('TAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTATTATATAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTAT'):
                print (line)
                f.write(f'{line}\n')        
         
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in mrep(Input[0]):
                print (line)
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
