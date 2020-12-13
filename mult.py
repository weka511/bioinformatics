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

#  MULT Multiple Alignment

import argparse
import os
import time
from   helpers import read_strings

from   align   import FindMultipleSequenceAlignment
 
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('MULT Multiple Alignment')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print(FindMultipleSequenceAlignment(['ATATCCG','TCCGA','ATGTACTG'],
                                                    score=lambda X: 1 if X[0]==X[1] and X[1]==X[2] and X[0]!='-' else 0))
        print (FindMultipleSequenceAlignment(['ATATCCG',
                             'TCCG',
                             'ATGTACTG',
                             'ATGTCTG']))
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        score, Alignment = FindMultipleSequenceAlignment([Input[1], Input[3], Input[5], Input[7]])
       
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            print (score)
            f.write(f'{score}\n')
            for line in Alignment:
                print (line)
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
