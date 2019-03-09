#    Copyright (C) 2019 Greenweaves Software Limited
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
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# GCON Global Alignment with Constant Gap Penalty

from align import align
from helpers import create_strings
from BA5J import san_kai

def get_indel_cost_with_skips(sigma,delta,i,j,di,dj,moves):
    i_previous = i+di
    j_previous = j+dj
    
    if (i_previous,j_previous) in moves:
        _,_,di_previous,dj_previous = moves[(i_previous,j_previous)]
        if (di_previous==0 and di==0) or (dj_previous==0 and dj==0):
            return i_previous-i,j_previous-j,delta

    return di,dj,sigma

def get_score(s,t,matrix,moves,showPath=False):
    return matrix[len(s)][len(t)],[],[]

def gcon(s,t):
    score,s1,t1 = san_kai([s0 for s0 in s],[t0 for t0 in t],sigma=5,epsilon=0)
    return score,''.join(s1),''.join(t1) 

if __name__=='__main__':
    from Bio.SubsMat.MatrixInfo import blosum62
    from Bio import SeqIO
    import sys,os
    if sys.argv[1]=='--sample':
        score,_,_=gcon('PLEASANTLY','MEANLY')
        #score,_,_=align('PLEASANTLY','MEANLY',
                            #replace_score=blosum62,
                            #indel_cost=(5,0),
                            #get_indel_cost=get_indel_cost_with_skips,
                            #backtrack=get_score,
                            #showScores=False,showPath=False)
        print (score)
    elif sys.argv[1]=='--test':
        inFile = open(r'C:\Users\Simon\Downloads\rosalind_gcon(3).txt','r')
        strings = []
        for record in SeqIO.parse(inFile,'fasta'):
            print (record.id)
            print (str(record.seq))
            strings.append(str(record.seq))        
        #strings=create_strings('gcon',fasta=True,ext=2)
        score,_,_=gcon(strings[0],strings[1])
        #score,_,_=align(strings[0],strings[1],
                        #replace_score=blosum62,
                        #indel_cost=(5,0),
                        #get_indel_cost=get_indel_cost_with_skips,
                        #backtrack=get_score,
                        #showScores=False,showPath=False)
        print (score)
    else:
        score,_,_=align(sys.argv[1],sys.argv[2],
                        replace_score=blosum62,
                        indel_cost=(5,0),
                        get_indel_cost=get_indel_cost_with_skips,
                        backtrack=get_score,
                        showScores=False,showPath=False)
        print (score)         
