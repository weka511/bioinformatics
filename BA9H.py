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

#  BA9H 	Pattern Matching with the Suffix Array 

import argparse
import os
import time
from   helpers import read_strings
from   snp import SuffixArray

def MatchOnePatternUsingSuffixArray(Text,Pattern,SuffixArray):
    minIndex = 0
    maxIndex = len(Text)
    while minIndex < maxIndex:
        midIndex = (minIndex + maxIndex)//2
        if Pattern>Text[SuffixArray[midIndex]:]:
            minIndex = midIndex + 1
        else:
            maxIndex = midIndex
            
    first    = minIndex
    maxIndex = len(Text)
    while minIndex<maxIndex:
        midIndex = (minIndex + maxIndex)//2
        if Text[SuffixArray[midIndex]:].startswith(Pattern):
            minIndex = midIndex + 1
        else:
            maxIndex = midIndex
            
    last = maxIndex
    return (first,last)

def MatchPatternsUsingSuffixArray(Text,Patterns):
    sfa = SuffixArray(Text)
    Matches = []
    for Pattern in Patterns:
        first,last = MatchOnePatternUsingSuffixArray(Text,Pattern,sfa)
        for i in range(first,last):
            Matches.append(sfa[i])
    return sorted(Matches)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9H 	Pattern Matching with the Suffix Array')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        #Text ='AATCGGGTTCAATCGGGGT'
        #sfa = SuffixArray(Text)
        #for i in range(len(sfa)):
            #print (i,sfa[i],Text[sfa[i]:])        
        #print (MatchOnePatternUsingSuffixArray(Text,'ATCG',sfa))
        #print (MatchOnePatternUsingSuffixArray(Text,'GGGT',sfa'ATCG'))
        #print (MatchOnePatternUsingSuffixArray(Text,'FOO',sfa))
        print (MatchPatternsUsingSuffixArray('AATCGGGTTCAATCGGGGT',['ATCG','GGGT']))
    

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = None
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
