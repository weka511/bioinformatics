#   Copyright (C) 2017-2020 Greenweaves Software Limited

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

# BINS Binary Search

import argparse
import os
import time
from helpers import read_strings

# bins
#
# Find a given set of keys in a given array.

def bins(n,m,a,values):
    
    def search(value,imin,imax):
        if imin>imax:
            return -1
        if imin==imax:
            return imin if value==a[imin-1] else -1
   
        imid = (imax+imin)//2
        left = search(value,imin,imid)
        
        return left if left>-1 else search(value,imid+1,imax)
            
    return [search(value,1,n) for value in values]


if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (bins(5,
                    6,
                    [10, 20, 30, 40, 50],
                    [40, 10, 35, 15, 40, 20]))
  
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Result = ' '.join([str(r) for r in bins(int(Input[0]),
                                                int(Input[1]),
                                                [int(i) for i in Input[2].split()],
                                                [int(i) for i in Input[3].split()])])
        
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s') 
    
