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
    
    def CreateProfile():
        Profile = [[0]*m for k in range(K)]
        for i in range(n):
            for j in range(m):
                ch    = Alignment[i][j]
                index = [k for k in range(len(Alphabet)) if Alphabet[k]==ch]
                if len(index)==1:
                    k = index[0]
                    Profile[k][j]+=1
                    
        for i in range(m):
            total = sum([Profile[k][i] for k in range(K)])
            for k in range(K):
                Profile[k][i] /= total
                 
        return Profile
    m          = len(Alignment[0])
    n          = len(Alignment)
    K          = len(Alphabet)
    Profile    = CreateProfile()
    States     = ['S','I0'] + flatten([CreateState(i) for i in range(m)]) + ['E']
    Transition = []
    Emission   = [[0]*m]
    return Transition, Emission

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA10E 	Construct a Profile HMM')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        Transition,Emission = ConstructProfileHMM(
                                   0.289,
                                   ['A',   'B',   'C',   'D',   'E'],
                                   ['EBA', 'EBD', 'EB-', 'EED', 'EBD', 'EBE','E-D','EBD'])
        
        print (Transition)
        print (Emission)
               
    if args.extra:
        Input,Expected                         = read_strings(f'data/ProfileHMM.txt',init=0)
        Transition,Emission = ConstructProfileHMM(
                                   float(Input[0]),
                                   Input[2].split(),
                                   Input[4:-1]) 
        print (Transition)
        print (Emission)
               
    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        Transition,Emission = ConstructProfileHMM(
                                   float(Input[0]),
                                   Input[2].split(),
                                   Input[4:-1]) 
        print (Transition)
        print (Emission)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for line in Result:
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
