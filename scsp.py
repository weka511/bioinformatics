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

# SCSP.py Interleaving Two Motifs

from align import align
from reference_tables import createSimpleDNASubst

def scsp(s,t):
      def extend(u1,u2,v1):
            w = []
            i = 0
            while u1[0]!=u2[i] and u1[0]!='-':
                  w.append(u2[i])
                  i+=1
            return w + u1,i*['-']+v1
      a,b,c = align(s,t,replace_score=createSimpleDNASubst(subst=len(s)+len(t)),indel_cost=0)
      b1,c1 = extend(b,s,c)
      c2,b2 = extend(c1,t,b1)
      super_sequence = [aa if bb=='-' else bb for aa,bb in zip(b2,c2)]
      return ''.join(super_sequence)      

if __name__=='__main__':
      from helpers import create_strings
      print (scsp('ATCTGAT','TGCATA'))
