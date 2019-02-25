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

def get_indel_cost(sigma,delta,i,j,di,dj,moves):
    i1   = i
    j1   = j
    di1  = di
    dj1  = dj
    cost = sigma
    while True:
        i1+=di
        j1+=dj
        if (i1,j1) in moves:
            _,_,di0,dj0 = moves[(i1,j1)]
            if (di0==0 and di1==0 and dj1!=0) or (dj0==0 and dj1==0 and di1!=0):
                di1+=di
                dj1+=dj
                continue

        return di1,dj1,cost

def get_score(s,t,matrix,moves,showPath=False):
    return matrix[len(s)][len(t)],[],[]

if __name__=='__main__':
    from Bio.SubsMat.MatrixInfo import blosum62
    from Bio import SeqIO
    score,_,_=align('PLEASANTLY','MEANLY',
                    replace_score=blosum62,indel_cost=(5,0),get_indel_cost=get_indel_cost,
                    backtrack=get_score,showScores=False,showPath=False)    
    #inFile = open(r'C:\Users\Simon\Downloads\rosalind_gcon.txt','r')
    #strings = []
    #for record in SeqIO.parse(inFile,'fasta'):
        #print (record.id)
        #print (str(record.seq))
        #strings.append(str(record.seq))    
    #score,_,_=align(strings[0],strings[1],
                    #replace_score=blosum62,indel_cost=(5,0),get_indel_cost=get_indel_cost,
                    #backtrack=get_score,showScores=False,showPath=False)
    print (score)
