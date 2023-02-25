#!/usr/bin/env python

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
# BA9F Find the Shortest Non-Shared Substring of Two Strings

import argparse
from helpers import read_strings
from snp import SuffixArray
import time
from numpy import argmin

def FindShortestNonShared(s,t):
    m = len(s)
    n = len(t)
    print (m,n)
    r0,p0,l0 = SuffixArray(s,auxiliary=True,padLCP=True)
    for i in range(len(r0)):
        print (f'{i:2d} {r0[i]:2d} {p0[i]:2d}  {l0[i]:2d} {s[r0[i]:]}')
    r1,p1,l1 = SuffixArray(t,auxiliary=True,padLCP=True)
    for i in range(len(r1)):
        print (f'{i:2d} {r1[i]:2d} {p1[i]:2d}  {l1[i]:2d} {t[r1q[i]:]}')
    #text = s + '$' + t + '#'
    #r,p,l = SuffixArray(text,auxiliary=True,padLCP=True)
    #Candidates  = []
    #Lengths     = []
    #previous    = -1
    #for i in range(len(r)):
        #if r[i]<m:
            #print (f'{i:2d} {r[i]:2d} {p[i]:2d}  {l[i]:2d} {text[r[i]:]}')
            #if i==previous+1:
                #Candidates.append(i)
                #Lengths.append(l[i])
            #previous = i
    #shortest = min(Lengths)
    #locations = [Candidates[i] for i in range(len(Candidates)) if Lengths[i]==shortest]
    #for i in locations:
        #print (f'{i:2d} {r[i]:2d} {p[i]:2d}  {l[i]:2d} {text[r[i]:]}')

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9D Find the Longest Repeat in a String')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (FindShortestNonShared('CCAAGCTGCTAGAGG','CATGCTGGGCTGGCT'))

    if args.extra:
        Input,Expected  = read_strings('data/ShortestNonSharedSubstring.txt',init=0)
        print (Expected[0])
        Actual          = FindShortestNonShared(Input[0],Input[1])
        print (len(Expected[0]),len(Actual))
        print (Expected[0])
        print (Actual)

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        Result = FindShortestNonShared(Input[0],Input[1])
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')

    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
