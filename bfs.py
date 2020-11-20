#    Copyright (C) 2017-2020 Greenweaves Software Limited
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
#    bfs 	Breadth First search

import argparse
import os
import time
from   helpers import read_strings
from   graphs  import ShortestDistances

def create_tree(links):
    max_node,_ = links[0]
    result     = {a:[] for a in range(1,max_node+1)}
    for (a,b) in links[1:]:
        result[a].append(b)
        
    return (max_node,result)

def format(Result):
    return ' '.join(str(i) for i in Result)
    
if __name__=='__main__':
    start  = time.time()
    parser = argparse.ArgumentParser('....')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print(
            format(
                ShortestDistances(
                    create_tree([(6, 6),
                                 (4, 6),
                                 (6, 5),
                                 (4, 3),
                                 (3, 5),
                                 (2, 1),
                                 (1, 4)]))))
    
    

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        Links  = [[int(x) for x in line.split()] for line in Input]
        Result = format(ShortestDistances(create_tree(Links)))
        print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            f.write(f'{Result}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
