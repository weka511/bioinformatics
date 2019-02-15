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

# LCSQ 	Finding a Shared Spliced Motif 

def lcsq_dp(s,t):
    lcs=[[0 for j in t] for i in s]
    def get_lcs(i,j):
        return lcs[i][j] if i>=0 and j>=0 else 0
    def populate():
        for i in range(len(s)):
            for j in range(len(t)):
                lcs[i][j] = (get_lcs(i-1,j-1) + 1) if s[i]==t[j] else max(get_lcs(i-1,j), get_lcs(i,j-1))
        for i in range(len(s)):
            print (lcs[i])
            
    def extract():        
        i = len(s)-1
        j = len(t)-1
        seq = []
        while True:
            if get_lcs(i,j-1) == get_lcs(i,j): j -= 1
            elif get_lcs(i-1,j) == get_lcs(i,j): i -= 1
            else:
                assert s[i]==t[j]
                seq.append(s[i])
                i -= 1
                j -= 1
            if  get_lcs(i,j)==0:
                print (i,j)
                return seq[::-1]
    populate()
    return extract()

def lcsq(s,t):
    def get_longest(seqs):
        best_so_far = []
        for seq in seqs:
            if len(seq)>len(best_so_far):
                best_so_far = seq
        return best_so_far
   
    def lcsq_helper(s,t):
        if len(s)==0: return t
        if len(t)==0: return s
        seqs=[]
        if (len(s)==1 or len(t)==1) and s[0]==t[0]:
            seqs.append([s[0]])
        for m in range(1,len(s)):
            for n in range(1,len(t)):
                seqs.append(lcsq_helper(s[0:m],t[0:n]) + lcsq_helper(s[m:],t[n:]))
        return get_longest(seqs)
    
    return ''.join(lcsq_dp([s0 for s0 in s], [t0 for t0 in t]))
    
if __name__=='__main__':
    print (lcsq('GTTCCTAATA','CGATAATTGAGA'))
    print (lcsq('AACCTTGG','ACACTGTGA'))
    
    print (lcsq_dp(['A','B','A','Z','D','C'],['B','A','C','B','A','D']))

