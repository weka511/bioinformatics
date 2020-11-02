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

#    NWC Negative Weight Cycle
#
#   Check whether a given graph contains a cycle of negative weight.
#   Given: A positive integer k, and k simple directed graphs with integer edge weights from -1000 to 1000
#   Return: For each graph, output 1 if it contains a negative weight cycle and -1 otherwise.

import argparse
import os
import time
from helpers import read_strings,create_list 
from graphs import bf  

def extract_data(Input):
    RowsAsInts = [[int(i) for i in row.split()] for row in Input]
    return [data if len(data)>1 else data[0] for data in RowsAsInts if len(data)>0]

          
def extract_graphs(data):
    Result = []
    Edges  = []
    for row in data[1:]:
        if len(row)==2:
            if len(Edges)>0:
                Result.append(Edges)
                Edges=[]
            Edges.append(row)
        else:
            assert len(row)==3
            Edges.append(row)
    Result.append(Edges)
    return Result

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('NWC Negative Weight Cycle')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        data = [2,
                 [4, 5],
                 [1, 4, 4],
                 [4, 2, 3],
                 [2, 3, 1],
                 [3, 1, 6],
                 [2, 1, -7],
                 
                 [3, 4],
                 [1, 2, -8],
                 [2, 3, 20],
                 [3, 1, -1],
                 [3, 2, -30]]
        
        print ([n for n,_,_ in [bf(edges) for edges in extract_graphs(data)]])

    if args.rosalind:
        Input     = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')      
        Result    = [n for n,_,_ in [bf(edges,s=0) for edges in extract_graphs(extract_data(Input))]] 
        Formatted = ' '.join(str(r) for r in Result)
        print (Formatted)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Formatted}\n')
                
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
