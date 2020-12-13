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

#  MULT Multiple Alignment

import argparse
import os
import time
from   helpers import read_strings
from   numpy   import argmax,ravel,zeros

def score(chars, match    = 0, mismatch = -1):    
    return sum(match if chars[i]==chars[j] else  mismatch for i in range(len(chars)) for j in range(i+1,len(chars)))

        
def FindHighestScoringMultipleSequenceAlignment(Strings,score = score):
    def indices(ns=[len(S)+1 for S in Strings]):
        N  = 1
        for n in ns:
            N *= n
        ii = [0]*len(ns)
        for _ in range(N-1):          
            for j in range(len(ns)-1,-1,-1):
                ii[j] += 1
                if ii[j]==ns[j]:
                    ii[j] = 0
                else:
                    break
            yield tuple(ii)
            
    def create_moves(m,options=[0,-1]):
        def create_raw_moves(m):
            if m==1:
                return [[o] for o in options]
            else:
                return [[o] + c for o in options for c in create_raw_moves(m-1)]
        return [move for move in tuple(create_raw_moves(m)) if len([x for x in move if x!=0])>0]
     
    def add(u,v):
        return tuple([a + b for (a,b) in zip(list(u),list(v))])
    

    def build_matrix():
        def calculate_scores(i):
            def get_score(move):
                previous = add(i,move)
                if len([p for p in previous if p<0]):
                    return None
                scorable = []
                for j in range(len(move)):
                    scorable.append(Strings[j][previous[j]]  if move[j]<0 else '-')
                return s[previous] + score(scorable)
            raw_scores = [(get_score(move),move) for move in available_moves]  # filter out -ve coordinates
            return [(r,m) for r,m in raw_scores if r != None]
        
        s     = zeros([len(S)+1 for S in Strings],dtype=int)
        path  = {}
        m     = len(Strings)
        available_moves = create_moves(m)
        for index_set in indices():
            scores_moves     = calculate_scores(index_set)
            scores           = [score for score,_ in scores_moves]
            moves            = [move  for _,move  in scores_moves]
            index_best_score = argmax(scores)
            s[index_set]     = scores[index_best_score]
            path[index_set]  = moves[index_best_score]
        return s,path,moves
    
    def backtrack(history):
        def reverse(S):
            return ''.join([s for s in S[::-1]])
        s,path,moves    = history
        position        = tuple([len(S) for S in Strings])
        alignment_score = s[position]
        Alignments      = [[] for S in Strings]
        while (len([p for p in position if p!=0]) >0):
            move = path[position]
            for j in range(len(move)):
                if move[j]==0:
                    Alignments[j].append('-')
                else:
                    Alignments[j].append(Strings[j][position[j]-1])
            position = add(position,move)
        return alignment_score,[reverse(s) for s in Alignments]
    return backtrack(build_matrix())





if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('MULT Multiple Alignment')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print(FindHighestScoringMultipleSequenceAlignment(['ATATCCG','TCCGA','ATGTACTG'],
                                                    score=lambda X: 1 if X[0]==X[1] and X[1]==X[2] and X[0]!='-' else 0))
        print (FindHighestScoringMultipleSequenceAlignment(['ATATCCG',
                             'TCCG',
                             'ATGTACTG',
                             'ATGTCTG']))
        
    

    if args.rosalind:
        Input  = read_strings(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')
 
        score, Alignment = FindHighestScoringMultipleSequenceAlignment([Input[1], Input[3], Input[5], Input[7]])
        #print (Result)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as f:
            print (score)
            f.write(f'{score}\n')
            for line in Alignment:
                print (line)
                f.write(f'{line}\n')
                
    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')    
