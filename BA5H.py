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
#    BA5H Find a Highest-Scoring Fitting Alignment of Two Strings

from reference_tables import createSimpleDNASubst
from rosalind_old import highest_scoring_global_alignment
from align import align
import numpy as np

def build_matrix(s,t,matrix,replace_score=createSimpleDNASubst(),indel_cost=1):
      moves = {}
      def score(pair):
            def reverse(pair):
                  a,b=pair
                  return (b,a)
            return replace_score[pair] if pair in replace_score else replace_score[reverse(pair)] 
      for i in range(len(s)+1):
            for j in range(len(t)+1):
                  if i==0 and j==0: next
                  elif i==0: 
                        matrix[i][j] = 0#matrix[i][j-1] - indel_cost
                        moves[(i,j)]  = (0,0-1,0,-1)
                  elif j==0:
                        matrix[i][j] = 0#matrix[i-1][j] - indel_cost
                        moves[(i,j)]  = (0,0,-1,0)
                  else:
                        scores       = [matrix[i-1][j]   - indel_cost,
                                matrix[i][j-1]   - indel_cost,
                                matrix[i-1][j-1] + score((s[i-1],t[j-1]))]
                        froms        = [(i-1,j,-1,0),
                                (i,j-1,0,-1),
                                (i-1,j-1,-1,-1)]
                        index        = np.argmax(scores)
                        matrix[i][j] = scores[index]
                        moves[(i,j)]  = froms[index]

      return matrix,moves      




def ba5h(s,t):
      d,s1,t1 =align([s0 for s0 in s], [t0 for t0 in t],build_matrix=build_matrix)
      return (d,''.join(s1),''.join(t1))

def ba5h_naive(s,t):
      best = -9999999
      s0   = ''
      t0   = ''
      for i in range(len(s)-len(t)):
            for j in range(len(t),2*len(t)):
                  score,s1,t1=highest_scoring_global_alignment(s[i:i+j],t,weights=createSimpleDNASubst(),sigma=1)
                  if score>best:
                        best = score
                        s0   = s1
                        t0   = t1
      return best,s0,t0

if __name__=='__main__':
      d,s1,t1 = ba5h('GTAGGCTTAAGGTTA','TAGATA')
      d,s1,t1 = ba5h('CAATCACCCCAATCCCTCAATCCTGGCCCCACGCATAGGCTAATGCCAATCGCGGCCAGGGTATAACCGCCATAACTGTGGGTCAGAAGGGATAAGTTCCACAATCCTATTTTCCTCGAGGCGCTTCGATGCGTTAACGCGTACACTCTGTCGGCCAACCGTGTGGGAGCCGAATTGGCTGGGCTGTTGAACATTCTATCAGTAGATAAACGAAGGTACATCCGAGGTTGTCGATCGACCGCGGGGTCGTAGCGCGTGCATGTTCCTTTCAGGCCCACATACTCCGGAACGGTTCATATCACGACTATTCTTGCACAATCGGACAACGGTGTACCATGGTGGACACCGTAGGAGACCAATACTGCGTAAATCATAAGCATTGGAGAGTGGACTGCTAGCGAGGCTCACCATGGAGTCTCGGTCGGCATCTCCTGACTGCTGTTCCATCGCGTTTTTCTTTTACTCACGCAATAAATCAATACCCCCTAACACAGGCCTGCTCCAGCCTTATTAAGGCCATAGTAGCTCTACATGTAGACCGAACGGAAGCACAGTTTGGTAGAAATTCTTAATCGACTATGGTCCGTGCAGGCCAAAAAAGGAATAATCTTCGAATTCTCACGCCTTCATTAGGGCGCACATGGTGGGGTAAATCACTGCACTCTGTTCGCAGTTAAGCGTTGCAATCAATATCGGCAGAACTCGGAGTCCGTATAAAGCCGCCTCAGCGTGCACACGCCCGTGCGGCACGTCATTAGACGAGGATTCCGGGGGACTGGCCTGTTCGTAATCCACTAAAACAATGGTCCTACCATCTAAAACGCACCGTGTTCCCCTCTACGGGAACCCCCTAGAT',
                     'AGAGCGCAGAGAAGTCATTAGAACATGTAGCACATCGCTTATTAAGGGTCAATACCTAAAGGGCCTAACTATACGCCACACGGAACAGCTC')     
      print (d)
      print (s1)
      print (t1)
# Expect:
# 22
# AGGGCGCACATG--GTGGGGTA-AATCAC-T-GCAC-TCTG-TTCGCAGTTAAGCGTTGCAATCAATATCGGC-AGAACTCGGAGTCCGTA--TAAAGCCGCCTCAGCGTGCACACGC-C
# AGAGCGCAGA-GAAGTCAT-TAGAA-CATGTAGCACATC-GCTT---A-TTAAG-G--G---TCAATA-C--CTA-AA---GG-G-CC-TAACTATA-C-GCCACA-CG-GAACA-GCTC
