# Copyright (C) 2019-2020 Greenweaves Software Limited

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

# SCC Strongly Connected Component

import argparse
import os
import time
from   helpers import read_strings, create_list
from   graphs  import scc 


if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('SCC Strongly Connected Component')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        nscc,adj,adj_R = scc([(6, 7),
                              (4, 1),
                              (1, 2),
                              (2, 4),
                              (5, 6),
                              (3, 2),
                              (5, 3),
                              (3, 5)
                              ])
        print (nscc,adj,adj_R)
        
   
    if args.rosalind:
        nscc,_,_ = scc(create_list(path='./data'))
        print (nscc)        
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{nscc}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
