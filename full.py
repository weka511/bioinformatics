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
    
    lefts = extract( [min(candidates.keys())],candidates)
    rights = extract([min([r for r in candidates.keys() if not r in lefts])],candidates)
    while len(lefts)>len(rights):
        for r in rights:
            if r in lefts:
                ii = lefts.index(r)
                lefts = lefts[:ii] + lefts[ii+1:]
    ll = [candidates[l][0][5] for l in lefts]
    return ''.join(ll[:-1])  

    
if __name__=='__main__':
    L = [          
        #1988.21104821, 
        #610.391039105,
        #738.485999105,
        #766.492149105,
        #863.544909105,
        #867.528589105,
        #992.587499105,
        #995.623549105,
        #1120.6824591,
        #1124.6661391,
        #1221.7188991,
        #1249.7250491,
        #1377.8200091
        14986.6971128,
        2717.10801638,
        2814.16077638,
        2880.17134638,
        2911.21353638,
        3009.21393638,
        3014.22272638,
        3145.26321638,
        3172.27726638,
        3232.29524638,
        3328.37837638,
        3360.35382638,
        3399.41548638,
        3507.42223638,
        3536.47439638,
        3610.43142638,
        3651.50133638,
        3723.51548638,
        3780.54392638,
        3820.56824638,
        3851.58103638,
        3967.63665638,
        3998.64944638,
        4080.72071638,
        4126.74440638,
        4195.74765638,
        4255.78699638,
        4326.78814638,
        4354.85540638,
        4425.85655638,
        4441.88743638,
        4538.94061638,
        4540.95584638,
        4639.98829638,
        4669.05080638,
        4726.07226638,
        4755.01523638,
        4854.13084638,
        4941.09454638,
        4982.22580638,
        5028.12657638,
        5085.23499638,
        5115.15860638,
        5199.27792638,
        5243.21718638,
        5362.34125638,
        5372.25977638,
        5465.35044638,
        5528.36088638,
        5593.40902638,
        5627.42929638,
        5708.43596638,
        5740.51335638,
        5795.46799638,
        5853.59741638,
        5940.62944638,
        5951.56910638,
        6011.66655638,
        6048.62186638,
        6112.71423638,
        6177.66445638,
        6213.76191638,
        6316.77110638,
        6363.74376638,
        6403.80313638,
        6460.82459638,
        6462.81217638,
        6533.84928638,
        6573.90865638,
        6719.92859638,
        6720.97706638,
        6777.99852638,
        6818.99700638,
        6849.03563638,
        6966.06541638,
        7035.11494638,
        7079.14947638,
        7136.16262638,
        7239.17181638,
        7265.22878638,
        7338.24022638,
        7379.27171638,
        7492.35577638,
        7494.34133638,
        7607.42539638,
        7648.45688638,
        7721.46832638,
        7747.52529638,
        7850.53448638,
        7907.54763638,
        7951.58216638,
        8020.63169638,
        8137.66147638,
        8167.70010638,
        8208.69858638,
        8265.72004638,
        8266.76851638,
        8412.78845638,
        8452.84782638,
        8523.88493638,
        8525.87251638,
        8582.89397638,
        8622.95334638,
        8669.92600638,
        8772.93519638,
        8809.03265638,
        8873.98287638,
        8938.07524638,
        8975.03055638,
        9035.12800638,
        9046.06766638,
        9133.09969638,
        9191.22911638,
        9246.18375638,
        9278.26114638,
        9359.26781638,
        9393.28808638,
        9458.33622638,
        9521.34666638,
        9614.43733638,
        9624.35585638,
        9743.47992638,
        9787.41918638,
        9871.53850638,
        9901.46211638,
        9958.57053638,
        10004.4713064,
        10045.6025664,
        10132.5662664,
        10231.6818764,
        10260.6248464,
        10317.6463064,
        10346.7088164,
        10445.7412664,
        10447.7564964,
        10544.8096764,
        10560.8405564,
        10631.8417064,
        10659.9089664,
        10730.9101164,
        10790.9494564,
        10859.9527064,
        10905.9763964,
        10988.0476664,
        11019.0604564,
        11135.1160764,
        11166.1288664,
        11206.1531864,
        11263.1816264,
        11335.1957764,
        11376.2656864,
        11450.2227164,
        11479.2748764,
        11587.2816264,
        11626.3432864,
        11658.3187364,
        11754.4018664,
        11814.4198464,
        11841.4338964,
        11972.4743864,
        11977.4831764,
        12075.4835764,
        12106.5257664,
        12172.5363364,
        12269.5890964
        
    ]
    
    print (full(L))  # expect KEKEP
