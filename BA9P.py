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

#  BA9P 	Implement TreeColoring 

import argparse
import os
import time
from   helpers import read_strings
from   snp     import ColourTree

            
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9P 	Implement TreeColoring ')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        Coloured = ColourTree({
                              0 : [],
                              1 : [],
                              2 : [0,1],
                              3 : [],
                              4 : [],
                              5 : [3,2],
                              6 : [],
                              7 : [4,5,6]},
                             {
                              0: 'red',
                              1: 'red',
                              3: 'blue',
                              4: 'blue',
                              6: 'red'})
        
        for node in sorted(Coloured.keys()):
            print (f'{node}: {Coloured[node]}')
    

    if args.rosalind:
        Input   = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        status  = 0
        adj     = {}
        colours = {}
        for line in Input:
            if line=='-':
                status = 1
                continue
            if status == 0:
                parts = line.split('->')
                adj[int(parts[0])] = [] if parts[1].strip()=='{}' else  [int(c) for c in parts[1].split(',')]               
            else:
                parts = line.split(':')
                colours[int(parts[0])] = parts[1].strip()               
            
        Result = ColourTree(adj,colours)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for key in sorted(Result.keys()):
                print (f'{key}: {Result[key]}')
                f.write(f'{key}: {Result[key]}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
