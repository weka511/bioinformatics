#  Copyright (C) 2020 Greenweaves Software Limited
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# BA9E Find the Longest Substring Shared by Two Strings

import argparse
import os
from   helpers import read_strings
from   snp     import SuffixArray
import sys
import time
from   numpy   import argmax

# LongestSharedSubstring
#
# Find the Longest Substring Shared by Two Strings
#
# See https://en.wikipedia.org/wiki/Longest_common_substring_problem

def LongestSharedSubstring(s,t):
            
    t0 = time.time()
    print ( f'About to build suffix Array {time.time()-t0}')
    text = s + '$' + t + '#'
    r,p,lcp = SuffixArray(text,auxiliary=True,padLCP=True)
    print ( f'Built tree {time.time()-t0}')
    
    previous_from_s = False
    Pairs = []
    for i in range(len(lcp)):
        from_s = r[i]<len(s)+2
        #print (f'{i},{r[i]},{lcp[i]},{text[r[i]:r[i]+lcp[i]]},{from_s},{previous_from_s}')
        if i>0 and previous_from_s != from_s:
            Pairs.append(i)
        previous_from_s = from_s
    candidate_LCPs = [lcp[i] if i in Pairs else 0 for i in range(len(lcp))]
    index = argmax(candidate_LCPs)
    return text[r[index]:r[index]+candidate_LCPs[index]]
        
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9E Find the Longest Substring Shared by Two Strings')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (LongestSharedSubstring('panama',
                                 'bananas'))
        print (LongestSharedSubstring('TCGGTAGATTGCGCCCACTC',
                                 'AGGGGCTCGCAGTGTAAGAA'))        
        
    if args.extra:
        Input,Expected  = read_strings('data/LongestSharedSubstring.txt',init=0)       
        Actual = LongestSharedSubstring(Input[0],Input[1])
        print (len(Expected[0]),len(Actual))
        print (Expected[0])
        print (Actual)        
 
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')    
        Result = LongestSharedSubstring(Input[0],Input[1])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')
            
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
