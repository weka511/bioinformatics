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

def lcsq(s,t):
    def get_prefix(s,t):
        prefix = []
        for i in range(min(len(s),len(t))):
            if s[i]==t[i]:
                prefix.append(s[i])
            else:
                break
        return prefix
    
    def lcsq_helper(s,t):
        if len(s)==0: return t
        if len(t)==0: return s
        seqs=[]
        for m in range(1,len(s)):
            for n in range(1,len(s)):
                seqs.append(lcsq_helper(s[0:m],t[0:n]) + lcsq_helper(s[m:],t[n:]))
        lcsq=[]
        for seq in seqs:
            if len(seq)>len(lcsq):
                lcsq=seq
        return lcsq
    
    s1 = [s0 for s0 in s]
    t1 = [t0 for t0 in t]
    prefix = get_prefix(s1,t1)
    suffix = list(reversed(get_prefix(list(reversed(s1)),
                                      list(reversed(t1)))))
    return ''.join(prefix + lcsq_helper(s1[len(prefix):len(s1)-len(suffix)],
                                        t1[len(prefix):len(t1)-len(suffix)]) + suffix)
    
if __name__=='__main__':
    print (lcsq('AACCTTGG','ACACTGTGA'))

