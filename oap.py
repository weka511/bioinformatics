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

#    OAP Overlap Alignment

import argparse
import os
import time
import sys
from   laff import read_fasta


def oap(v,w,
        match_bonus   = 1,
        mismatch_cost = 2,
        indel_cost    = 2):
    def score(v_,w_):
        return match_bonus if v_==w_ else -mismatch_cost
        
    def dynamic_programming(v,w):
        def format(i,j):
            return v[i-1] if j==0 else f'{distances[i][j]}'
        distances = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
        for i in range(1,len(v)+1):
            for j in range(1,len(w)+1):
                distances[i][j]  = max(
                                        distances[i-1][j]   - indel_cost,
                                        distances[i][j-1]   - indel_cost,
                                        distances[i-1][j-1] + score(v[i-1],w[j-1]))
        
        #print ('  ' + ' '.join(w[j] for j in range(len(w))))        
        #for i in range(1,len(v)+1): 
            #print (' '.join(format(i,j) for j in range(len(w)+1)))
        i = len(v)
        distance = max(distances[i])
        for j in range(len(w)+1):
            if distances[i][j] == distance: break
        #distance,i,j = get_max_score(distances, len(v), len(w))
        print (distance,i,j)
        v1       = []
        w1       = []
        while i>0 and j>0:
            #if distances[i][j]==0:
                #i1,j1 = (0,j-1)
            if distances[i][j]==distances[i-1][j]   - indel_cost:
                i1,j1 = (i-1,j)
            elif distances[i][j]==distances[i][j-1]   - indel_cost:
                i1,j1 = (i,j-1)
            elif distances[i][j]==distances[i-1][j-1] + score(v[i-1],w[j-1]):
                i1,j1 = (i-1,j-1)
            else:
                raise Exception(f'This cannot possible happen {i} {j}!')
             
            v1.append(v[i1] if i1<i else '-')
            w1.append(w[j1] if j1<j else '-')
            #if distances[i][j]==0: break
            i,j = i1,j1
    
        return distance,v1[::-1],w1[::-1] # was sum(score(u,v) for u,v in zip(v1,w1))
    
    score,u1,v1=dynamic_programming([vv for vv in v],[ww for ww in w])
    return score,''.join(u1),''.join(v1)

if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('OAP Overlap Alignment')
    parser.add_argument('--sample',    default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--rosalind',  default=False, action='store_true', help='process Rosalind dataset')
    parser.add_argument('--version',   default=False, action='store_true', help='Get version of python')
    parser.add_argument('--frequency', default=100,   type=int,            help='Number of iteration per progress tick' )
    args = parser.parse_args()
    
    if args.version:
        print (f'{sys.version}')
        
    if args.sample:
        d,s1,t1 = oap('CTAAGGGATTCCGGTAATTAGACAG',
                      'ATAGACCATATGTCAGTGACTGTGTAA')
        
        print ('{0}'.format(d))
        print (s1)
        print (t1)
        
    if args.rosalind:
        Data      = read_fasta(f'data/rosalind_{os.path.basename(__file__).split(".")[0]}.txt')          
        d,s1,t1   = oap(Data[0],Data[1])      
        print ('{0}'.format(d))
        print (s1)
        print (t1)
        with open(f'{os.path.basename(__file__).split(".")[0]}.txt','w') as o:
            o.write('{0}\n'.format(d))
            o.write('{0}\n'.format(s1))
            o.write('{0}\n'.format(t1))    

    elapsed = time.time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')        