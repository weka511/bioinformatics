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

from helpers import create_strings
from BA5J import san_kai

def gcon(s,t):
    score,s1,t1 = san_kai([s0 for s0 in s],[t0 for t0 in t],sigma=5,epsilon=0)
    return score,''.join(s1),''.join(t1) 

if __name__=='__main__':
    from Bio.SubsMat.MatrixInfo import blosum62
    import sys
    score = -float('inf')
    if sys.argv[1]=='--sample':
        score,_,_=gcon('PLEASANTLY','MEANLY')
    elif sys.argv[1]=='--test':
        strings   = create_strings('gcon',ext=3,fasta=True)    
        score,_,_ = gcon(strings[0],strings[1])
    else:
        score,_,_=gcon(sys.argv[1],sys.argv[2])
    print (score)         
