#!/usr/bin/env python

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

#  BA10H 	Estimate the Parameters of an HMM

import argparse
import os
import time
from   helpers import read_strings
from   hmm     import float2str,EstimateParameters



# formatEmission

def formatEmission(Emission,States,Alphabet,precision=2):
    yield '\t'+'\t'.join(Alphabet)
    for state in States:
        row = []
        for symbol in Alphabet:
            probability = Emission[(symbol,state)]
            row.append(float2str(probability,precision))
        yield state + '\t' + '\t'.join(row)


# formatTransition



def formatTransition(Transition,States,precision=2):
    yield '\t'+ '\t'.join(state for state in States)
    for state1 in States:
        row = []
        for state2 in States:
            probability = Transition[(state1,state2)] if (state1,state2) in Transition else 0
            row.append(float2str(probability,precision))
        yield state1 + '\t' + '\t'.join(row)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA10H 	Estimate the Parameters of an HMM ')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',     default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--precision', default=3,                          help='Controls display of probabilities')
    args = parser.parse_args()
    if args.sample:
        Transitions,Emissions = EstimateParameters('yzzzyxzxxx',
                           ['x',  'y',   'z'],
                           'BBABABABAB',
                           ['A', 'B', 'C'])
        for row in formatTransition(Transitions,['A', 'B', 'C'],precision=args.precision):
            print (row)
        for row in formatEmission(Emissions,['A', 'B', 'C'], ['x',  'y',   'z'],precision=args.precision):
            print (row)

    if args.extra:
        Input,Expected             = read_strings(f'data/HMMParameterEstimation.txt',init=0)
        Transitions,Emissions = EstimateParameters(Input[0],Input[2].split(),Input[4],Input[6].split())
        for row in formatTransition(Transitions,Input[6].split(),precision=args.precision):
            print (row)
        print ('--------')
        for row in formatEmission(Emissions,Input[6].split(), Input[2].split(),precision=args.precision):
            print (row)

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
        Transitions,Emissions = EstimateParameters(Input[0],Input[2].split(),Input[4],Input[6].split())

        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            for row in formatTransition(Transitions,Input[6].split(),precision=args.precision):
                print (row)
                f.write(f'{row}\n')
            print ('--------')
            f.write('--------\n')
            for row in formatEmission(Emissions,Input[6].split(), Input[2].split(),precision=args.precision):
                print (row)
                f.write(f'{row}\n')

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
