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

#  BA10E Construct a Profile HMM

import abc
import argparse
import os
import time
from   helpers import read_strings
from   hmm import ConstructProfileHMM,formatEmission,formatTransition


        
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA10E Construct a Profile HMM')
    parser.add_argument('--sample',    default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',     default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind',  default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--text',      default=False, action='store_true', help='process dataset from textbook')
    parser.add_argument('--precision', default=3,                          help='Controls display of probabilities')

    args = parser.parse_args()
    if args.sample:
        States,Transition,Emission = ConstructProfileHMM(0.289,
                                                         ['A',   'B',   'C',   'D',   'E'],
                                                         ['EBA', 'EBD', 'EB-', 'EED', 'EBD', 'EBE','E-D','EBD'])
        
        for row in formatTransition(Transition,States,precision=args.precision):
            print (row)
        print ('--------')
        for row in formatEmission(Emission,States,['A',   'B',   'C',   'D',   'E'],precision=args.precision):
            print (row)
            
    if args.text:
        Alphabet = ['A',   'C',   'D',   'E', 'F']
        States,Transition,Emission = ConstructProfileHMM(0.35,
                                                         Alphabet,
                                                         ['ACDEFACADF',
                                                          'AFDA---CCF',
                                                          'A--EFD-FDC',
                                                          'ACAEF--A-C',
                                                          'ADDEFAAADF'])
        
        for row in formatTransition(Transition,States,precision=args.precision):
            print (row)
        print ('--------')
        for row in formatEmission(Emission,States,Alphabet,precision=args.precision):
            print (row)
            
    def compare(Expected,Actual):
        Mismatch = False
        if len(Expected)!=len(Actual):
            print(f'{len(Expected)}!={len(Actual)}')
            Mismatch = True
        for i in range(min(len(Expected),len(Actual))):
            if Expected[i] != Actual[i]:
                print (f'{i} {Expected[i]} {Actual[i]}')
                Mismatch = True
        return Mismatch
    
    if args.extra:
        Input,Expected             = read_strings(f'data/ProfileHMM.txt',init=0)
        States,Transition,Emission = ConstructProfileHMM(float(Input[0]),
                                                         Input[2].split(),
                                                         Input[4:-1])
        i = 0
        mismatches = 0
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for row in formatTransition(Transition,States,precision=args.precision):
                print (row)
                f.write(f'{row}\n')
                if compare(Expected[i],row):
                    mismatches+=1
                i += 1
            print ('--------')
            f.write('--------\n') 
            if compare(Expected[i],'--------'):
                mismatches+=1
            i += 1            
            for row in formatEmission(Emission,States,Input[2].split(),precision=args.precision):
                print (row)
                f.write(f'{row}\n')
                if compare(Expected[i],row):
                    mismatches+=1
                i += 1
        print (mismatches)
               
    if args.rosalind:
        Input                      = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        States,Transition,Emission = ConstructProfileHMM(float(Input[0]),
                                                         Input[2].split(),
                                                         Input[4:-1]) 

        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for row in formatTransition(Transition,States,precision=args.precision):
                print (row)
                f.write(f'{row}\n')
            print ('--------')
            f.write('--------\n')
            for row in formatEmission(Emission,States,Input[2].split(),precision=args.precision):
                print (row)            
                f.write(f'{row}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
