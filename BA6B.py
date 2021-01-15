#    Copyright (C) 2019-2021 Greenweaves Software Limited
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>
#
#    BA6B Compute the Number of Breakpoints in a Permutation 

import argparse
import os
import time
from   helpers import read_strings
from fragile import getBreakPoints,get_permutation



if __name__=='__main__':
     start = time.time()
     parser = argparse.ArgumentParser('Find a Shortest Transformation of One Genome into Another by 2-Breaks')
     parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
     parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
     args = parser.parse_args()
     if args.sample:   
          print (getBreakPoints([+3, +4, +5, -12, -8, -7, -6, +1, +2, +10, +9, -11, +13, +14]))
          
     if args.rosalind:
          Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
          Result = getBreakPoints(get_permutation(Input[0]))
          print (Result)
          
     elapsed = time.time() - start
     minutes = int(elapsed/60)
     seconds = elapsed - 60*minutes
     print (f'Elapsed Time {minutes} m {seconds:.2f} s')       