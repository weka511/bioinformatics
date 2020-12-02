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

#  BA9M 	Implement BetterBWMatching 

import argparse
import os
import time
from   helpers import read_strings
#from   snp     import BW_Match


def BetterBWMatching(LastColumn,Patterns):
    
    def Count(symbol,i):
        return sum([1 for ch in LastColumn[:i] if ch==symbol]) 

    def FirstOccurence(ch):
        for i in range(len(FirstColumn)):
            if FirstColumn[i]==ch:
                return i
     
    def LastColumnContains(symbol,top,bottom):
        topIndex    = None
        bottomIndex = None
        N           = len(LastColumn)
        
        for i in range(top,N):
            if LastColumn[i]==symbol:
                bottomIndex = i
                break
            
        for i in range(bottom,-1,-1):
            if LastColumn[i]==symbol:
                topIndex = i
                break
            
        return topIndex,bottomIndex
    
    def Match(Pattern):
        top    = 0
        bottom = len(LastColumn) - 1
        while top <= bottom:
            if len(Pattern) > 0:
                symbol  = Pattern[-1]
                Pattern = Pattern[:-1]
                topIndex,bottomIndex = LastColumnContains(symbol,top,bottom)
                if type(topIndex)==int and type(bottomIndex)==int:
                    top    = FirstOccurences[symbol] + Count(symbol,top)
                    bottom = FirstOccurences[symbol] + Count(symbol,bottom+1) - 1
                else:
                    return 0
            else:
                return bottom - top + 1
        return 0    
    
 
    FirstColumn = sorted(LastColumn)
    FirstOccurences = {ch:FirstOccurence(ch) for ch in FirstColumn}
    return [Match(Pattern) for Pattern in Patterns]

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9M 	Implement BetterBWMatching ')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (BetterBWMatching('GGCGCCGC$TAGTCACACACGCCGTA',
                                ['ACC', 'CCG', 'CAG']))
        
    

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
