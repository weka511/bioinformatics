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
#    full Inferring Peptide from Full Spectrum

from spectrum import create_lookup,get_abbrev
from reference_tables import amino_acids
from rosalind import get_weight

def get_n(L):
    n = (len(L)-3)//2
    assert(2*n+3==len(L))
    return n

#    full
#
#    Inferring Peptide from Full Spectrum
#
#    Inputs : A list L containing 2n+3 positive real numbers (n<=100). 
#             The first number in L is the parent mass of a peptide P, and all 
#             other numbers represent the masses of some b-ions and y-ions of P
#             (in no particular order). You may assume that if the mass of a b-ion is present, 
#             then so is that of its complementary y-ion, and vice-versa.
#
#    Return: A protein string t
#            of length n for which there exist two positive real numbers w1 and w2
#            such that for every prefix p and suffix s of t, each of w(p)+w1 and w(s)+w2
#            is equal to an element of L.
#            (In other words, there exists a protein string whose t-prefix and 
#            t-suffix weights correspond to the non-parent mass values of L.)
#            If multiple solutions exist, you may output any one.

def full(s,epsilon=0.000001):
    def extract(seq,candidates):
        while True:
            key = seq[-1]
            if not key in candidates:
                return seq
            succs = candidates[key]
            i,j,l1,l2,diff,candidate,_ = succs[0]
            seq.append(j)    
    masses,pairs   = create_lookup()
    n              = get_n(L)
    diffs          = [(i,j,L[i],L[j],abs(L[i] - L[j])) for i in range(1,len(L)) for j in range(i+1,len(L))]
    candidates     = {}
    for i,j,l1,l2,diff in diffs:
        abbrev    = get_abbrev(diff,masses,pairs)
        candidate = abbrev if abs(diff-amino_acids[abbrev].mon_mass)<epsilon else None
        if candidate != None:
            if not i in candidates:
                candidates[i]=[]
            candidates[i].append((i,j,l1,l2,diff,candidate,abs(diff-amino_acids[abbrev].mon_mass)))
    for key in candidates:
        candidates[key]= sorted(candidates[key],key=lambda x:x[6])
  #      print (key,candidates[key])
    
    lefts = extract( [min(candidates.keys())],candidates)
    rights = extract([min([r for r in candidates.keys() if not r in lefts])],candidates)
    while len(lefts)>len(rights):
        for r in rights:
            if r in lefts:
                ii = lefts.index(r)
                lefts = lefts[:ii] + lefts[ii+1:]
    ll = [candidates[l][0][5] for l in lefts]
    return ''.join(ll[:-1])  
##   create_candidates
##
##   Find list of amino acids that could possibly explain some of spectrum
    #def create_candidates():
        #masses,pairs = create_lookup()
        #candidates=[]
        #for diff in [m1-m0 for m0 in s[1:] for m1 in s[1:] if m1>m0]:
            #a = get_abbrev(diff,masses,pairs)
            #m = amino_acids[a].mon_mass
            #if abs(diff-amino_acids[a].mon_mass)<epsilon:
                #candidates.append(a)
        #return candidates
    
    ## create_pairs
    ##
    ## Extract a list of pairs (x,y) such that x+y == s[0].
    ## Each pair is either of form (b_ion,y_ion) or (y_ion,b_ion)
    
    #def create_pairs():
        #return [(s[i],s[-i]) for i in range(1,len(s)//2+1)]

    #def lookup(diff,fmasses,pairs):
        #a = get_abbrev(diff,masses,pairs)
        #m = amino_acids[a].mon_mass
        #if abs(diff-amino_acids[a].mon_mass)<epsilon:
            #return a
    
    #masses,stuff = create_lookup()
    #pairs = create_pairs()
    #print (len(pairs),pairs)
    #for (a,b) in pairs[1:]:
        #x=lookup(a-pairs[0][0],masses,stuff)
        #y=lookup(b-pairs[0][0],masses,stuff)
        #print (a-pairs[0][0],b-pairs[0][0],x,y)
    ##candidates= create_candidates()
    ##print (len(candidates),candidates)
    #protein = []
    #return ''.join(protein)

#def lluf(peptide,L):
    #masses,pairs   = create_lookup()
    #n              = get_n(L)
    #suffixes       = [peptide[i:] for i in range(len(peptide))]
    #prefixes       = [peptide[:i] for i in range(1,len(peptide))]
    #suffix_weights = [get_weight(s) for s in suffixes]
    #prefix_weights = [get_weight(p) for p in prefixes]
    #diffs = [(i,j,L[i],L[j],abs(L[i] - L[j])) for i in range(1,len(L)) for j in range(i+1,len(L))]
    #for i,j,l1,l2,d in diffs:
        #print (i,j,l1,l2,d)
    #return []
    #ll=[get_weight(s[i:]) for i in range(len(s))]
    #rr=[get_weight(s[:i]) for i in range(1,len(s))]
    #return sorted(ll+rr)
    
if __name__=='__main__':
    L = [          
        1988.21104821, 
        610.391039105,
        738.485999105,
        766.492149105,
        863.544909105,
        867.528589105,
        992.587499105,
        995.623549105,
        1120.6824591,
        1124.6661391,
        1221.7188991,
        1249.7250491,
        1377.8200091
    ]
    
    #print (lluf('KEKEP',L))
    print (full(L))  # expect KEKEP
