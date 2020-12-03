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

#  BA10E 	Construct a Profile HMM

import argparse
import os
import time
from   helpers import read_strings, flatten

def ConstructProfileHMM(theta,Alphabet,Alignment):
    def CreateState(i):
        return [f'{c}{i+1}' for c in 'MDI']
    
    def CreateAllStates():
        return ['S','I0'] + flatten([CreateState(i) for i in range(m)]) + ['E']
    
    m          = len(Alignment[0])
    States     = CreateAllStates()
    Transition = []
    Emission   = []
    x=0

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA10E 	Construct a Profile HMM')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        ConstructProfileHMM(0.289,
                            ['A',   'B',   'C',   'D',   'E'],
                            ['EBA', 'EBD', 'EB-', 'EED', 'EBD', 'EBE','E-D','EBD'])
        
    
    if args.extra:
        Input,Expected                         = read_strings(f'data/ProfileHMM.txt',init=0)
        ConstructProfileHMM(float(Input[0]), Input[2].split(),Input[4:-1]) 
        
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
