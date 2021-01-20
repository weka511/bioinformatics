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
#    EDTA Edit Distance Alignment http://rosalind.info/problems/edta/


import argparse
import os
import time
from   helpers import read_strings
from   Bio import SeqIO
from   align import edta

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('EDTA Edit Distance Alignment')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    
    if args.sample:
        d,s,t=edta('PRETTY' ,'PRTTEIN')
        print(d)
        print (s)
        print (t)
        
    if args.rosalind:
        inFile = open(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt','r')
        strings = []
        for record in SeqIO.parse(inFile,'fasta'):
            strings.append(str(record.seq))
        d,s,t=edta(strings[0],strings[1])
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            print(d)
            f.write(f'{d}\n')
            print (s)
            f.write(f'{s}\n')
            print (t)
            f.write(f'{t}\n')
        
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')        
