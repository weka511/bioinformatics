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
#
#    prsm Matching a Spectrum to a Protein

from spectrum import conv
from reference_tables import amino_acids
from numpy import argmax

def complete_spectrum(P):
    def spectrum(S):
        return sum([amino_acids[s].mon_mass for s in S])
    prefixes = [P[:i] for i in range(1,len(P))] +[P]
    suffixes = [P[i:] for i in range(1,len(P))]
    ss= [spectrum(p) for p in prefixes + suffixes]
    return ss#[item for sublist in ss for item in sublist]



def prsm(s,R):
    def count(c):
        return c[1][0]
    m = 0
    i = 0
    Ss = [(P,complete_spectrum(P)) for P in s]
    #for p,ss in Ss:
        #print (p,ss)
    Cs = [(P,conv(R,S1,eps=0.00001)) for (P,S1) in Ss]
    #for (p,(a,b)) in Cs:
        #print (p,a,b)
    ii = argmax([count(c) for c in Cs])
    #print (ii,Cs[ii])
    a,(b,c) = Cs[ii]
    return b,s[ii]

if __name__=='__main__':
    from helpers import create_strings
    #m,s_max = prsm(['GSDMQS',
                    #'VWICN',
                    #'IASWMQS',
                    #'PVSMGAD'],
                   #[445.17838,
                    #115.02694,
                    #186.07931,
                    #314.13789,
                    #317.1198,
                    #215.09061]
                   #)

    
    i = 0
    s = []
    R = []
    for ll in create_strings(ext=3):
        if i ==0:
            n = int(ll)
        elif i<n+1:
            s.append(ll)
        else:
            R.append(float(ll))
        i+=1
            
    m,s_max = prsm(s,R)
    
    print (m)
    print (s_max)    
    