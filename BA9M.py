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
from   snp     import BetterBWMatching


if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9M 	Implement BetterBWMatching ')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (BetterBWMatching('GGCGCCGC$TAGTCACACACGCCGTA',
                                ['ACC', 'CCG', 'CAG']))
        
    
    if args.extra:
        Input,Expected  = read_strings(f'data/BetterBWMatching.txt',init=0)
        Expected = [int(s) for s in Expected[0].split()]
        Result   = BetterBWMatching(Input[0],Input[1].split())
        if len(Result)!=len(Expected):
            print (f'Mismatched lengths {len(Result)} {len(Expected)}')
            
        mismatches = 0
        for i in range(min(len(Result),len(Expected))):
            if (Result[i]!=Expected[i]):
                print (f'{i} {Expected[i]} {Result[i]}')
                mismatches += 1
        print (f'{mismatches} mismatches')
            
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = BetterBWMatching(Input[0],Input[1].split())
        print (Result)
        line = ' '.join(str(i) for i in Result)
        print (line)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
