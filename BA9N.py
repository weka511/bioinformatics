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

# BA9N 	Find All Occurrences of a Collection of Patterns in a String 

import argparse
import os
import time
from   helpers import read_strings
from   snp import getFirstOccurrence, getCount, ColumnContains, BurrowsWheeler,PartialSuffixArray,LastToFirst,getN,get_char

def EvenBetterBWMatching(Text,Patterns,K=10):
    def find_position(i):
        steps = 0
        while i%K != 0:
            ch = LastColumn[i]
            n   = getN(ch,i,LastColumn)  
            i   = get_char(ch,FirstColumn,n)
            steps += 1
        return steps + i
    
    # Match
    #
    # Used by BetterBWMatching to match one occurrence of Pattern in text 
    #
    # Parameters:
    #      Pattern
    #      FirstOccurrences
    def Match(Pattern,FirstOccurrences):
        top    = 0
        bottom = len(LastColumn) - 1
        while top <= bottom:
            if len(Pattern) > 0:
                symbol  = Pattern[-1]
                Pattern = Pattern[:-1]
                topIndex,bottomIndex = ColumnContains(LastColumn,symbol,top,bottom)
                if type(topIndex)==int and type(bottomIndex)==int:
                    top    = FirstOccurrences[symbol] + getCount(LastColumn,symbol,top)
                    bottom = FirstOccurrences[symbol] + getCount(LastColumn,symbol,bottom+1) - 1
                else:
                    return []
            else:
                return [find_position(i) for i in range(top,bottom+1)]
 
        return []   
    
    PSA              = [i for _,i in PartialSuffixArray(Text+'$',K)]
    LastColumn       = BurrowsWheeler(Text+'$')
    FirstColumn      = sorted(LastColumn)
    FirstOccurrences = {ch:getFirstOccurrence(FirstColumn,ch) for ch in FirstColumn}
    return [match for Pattern in Patterns for match in Match(Pattern,FirstOccurrences)]

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9N Find All Occurrences of a Collection of Patterns in a String ')
    parser.add_argument('--panama',   default=False, action='store_true', help='process panama bananas')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',   default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    
    if args.panama:
        print (EvenBetterBWMatching('panamabananas', ['ana'], K=5))  
        
    if args.sample:
        print (EvenBetterBWMatching('AATCGGGTTCAATCGGGGT', ['ATCG', 'GGGT'], K=5))
        
    if args.extra:
        Input,Expected  = read_strings('data/MultiplePatternMatching.txt',init=0)
        Actual = EvenBetterBWMatching(Input[0],Input[1:])
        print (sorted(Expected[0]))
        print (sorted(Actual))
        
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in EvenBetterBWMatching(Input[0],Input[1:]):
                print(line)
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
