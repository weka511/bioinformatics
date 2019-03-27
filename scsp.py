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

# scsp
#
# Find a shortest common supersequence to two sequences
#
# Paramaters: a One sequence as a string
#             b The other sequence
#
# Returns: String - shortest common supersequence to a and b
#
# Strategy: find an alignment, using a scoring matrix that penalizes mismatches severely,
#           but does not penalize indels. Then extend because there may be a few missed characters
#           at the start of either string

def scsp(s,t):
      #  extend
      #
      #  Add missing characters to start
      #
      #  Parameters: extend_me  List from alignment
      #              original   Original before alignmnet
      #              pad_me     Companion - will be padded with leading spaces
      #
      # Returns: two strings, the first extended, the second padded
      
      def extend(extend_me,original,pad_me):
            extension = []
            i = 0
            while extend_me[0]!=original[i] and extend_me[0]!='-':
                  w.append(original[i])
                  i+=1
            return extension + extend_me,i*['-']+pad_me
      
      _,b,c          = align(s,t,replace_score=createSimpleDNASubst(subst=len(s)+len(t)),indel_cost=0)
      b1,c1          = extend(b,s,c)
      c2,b2          = extend(c1,t,b1)
      super_sequence = [aa if bb=='-' else bb for aa,bb in zip(b2,c2)]
      return ''.join(super_sequence)      

if __name__=='__main__':
      from helpers import create_strings
      print (scsp('TCCAGCGTCCTTAGGGCCGGTGCCGGTACTCCTAGGGAGATTTCTTGGGACCAATGCCACAGACAAGGTTGCGACAACCATACGGGCCGAAACCTA',
                  'TTCTTAGGGTCAGTCTGAAATGCGGGAACCATCCGTTATATCATCACGTCTGTCGTATTGGCCTGGCCCCTACGGTTGGCTG'))
