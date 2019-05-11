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

# OAP Overlap Alignment

from numpy import zeros,argmax,dtype,int,shape,unravel_index
from time import time
from reference_tables import createSimpleDNASubst

def oap(v,w,match_bonus=+1,mismatch_cost=2,indel_cost=2):
    def dynamic_programming(v,w):
        distances = zeros((len(v)+1,len(w)+1),dtype=int)
 
        for i in range(1,len(v)+1):
            for j in range(1,len(w)+1):
                distances[i,j]  = max( distances[i-1,j]   - indel_cost,
                                       distances[i,j-1]   - indel_cost,
                                       distances[i-1,j-1] + (match_bonus if v[i-1]==w[j-1] else -mismatch_cost))
                
        i,j     = unravel_index(distances.argmax(), distances.shape)
        distance = distances[i,j]
        v1       = []
        w1       = []
        while i>0 and j>0:
            index =  argmax([distances[i-1,j]   - indel_cost,
                             distances[i,j-1]   - indel_cost,
                             distances[i-1,j-1] + (match_bonus if v[i-1]==w[j-1] else -mismatch_cost)])
            i1,j1 = [(i-1,j),(i,j-1),(i-1,j-1)][index]
            v1.append(v[i1] if i1<i else '-')
            w1.append(w[j1] if j1<j else '-')
            i,j=i1,j1
    
        return distance,v1[::-1],w1[::-1]
    
    score,u1,v1=dynamic_programming([vv for vv in v],[ww for ww in w])
    return score,''.join(u1),''.join(v1)

if __name__=='__main__':
    from helpers import create_strings
    #d,s1,t1 = oap('CTAAGGGATTCCGGTAATTAGACAG','ATAGACCATATGTCAGTGACTGTGTAA')
    start_time = time()
    strings  = create_strings('oap',fasta=1,ext=3)
    d,s1,t1 = oap(strings[0],strings[1])      
    print ('{0}'.format(d))
    print (s1)
    print (t1)
    print (time()-start_time)
    with open('oap.txt','w') as o:
        o.write('{0}\n'.format(d))
        o.write('{0}\n'.format(s1))
        o.write('{0}\n'.format(t1))    
    