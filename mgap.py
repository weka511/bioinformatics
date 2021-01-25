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

#  MGAP Maximizing the Gap Symbols of an Optimal Alignment

import argparse
import os
import time
from   helpers import read_strings

# mgap
#
# For the computation of an alignment score generalizing the edit alignment score, 
# let m denote the score assigned to matched symbols, d denote the score assigned to mismatched non-gap symbols,
# and g denote the score assigned a symbol matched to a gap symbol '-' (i.e., gis a linear gap penalty).

# Given: Two DNA strings s and t,  each of length at most 5000 bp.
#
# Return: The maximum number of gap symbols that can appear in any maximum score alignment of s
# and t with score parameters satisfying m>0, d<0, and g<0.

def mgap(s,t):
    pass

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('MGAP Maximizing the Gap Symbols of an Optimal Alignment')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (mgap('AACGTA','ACACCTA'))
        
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
