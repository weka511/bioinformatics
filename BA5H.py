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

import numpy as np
from rosalind_old import highest_scoring_global_alignment



def ba5h(s,t,show_scores=False):
      def createWeights():
            bases=['A', 'T', 'G', 'C']
            weights={}
            for i in range(len(bases)):
                  for j in range(len(bases)):
                        weights[(bases[i],bases[j])]=+1 if i==j else -1           
            return weights
      best = -9999999
      s0   = ''
      t0   = ''
      
      for i in range(len(s)):
            for j in range(len(t),2*len(t)):
                  #print (i,j)
                  score,s1,t1=highest_scoring_global_alignment(s[i:j],t,weights=createWeights(),sigma=1)
                  if score>best:
                        print (score)
                        best = score
                        s0   = s1
                        t0   = t1
      return best,s0,t0

if __name__=='__main__':
#      d,s1,t1 = ba5h('GTAGGCTTAAGGTTA','TAGATA',show_scores=True)
      d,s1,t1 = ba5h('CAATCACCCCAATCCCTCAATCCTGGCCCCACGCATAGGCTAATGCCAATCGCGGCCAGGGTATAACCGCCATAACTGTGGGTCAGAAGGGATAAGTTCCACAATCCTATTTTCCTCGAGGCGCTTCGATGCGTTAACGCGTACACTCTGTCGGCCAACCGTGTGGGAGCCGAATTGGCTGGGCTGTTGAACATTCTATCAGTAGATAAACGAAGGTACATCCGAGGTTGTCGATCGACCGCGGGGTCGTAGCGCGTGCATGTTCCTTTCAGGCCCACATACTCCGGAACGGTTCATATCACGACTATTCTTGCACAATCGGACAACGGTGTACCATGGTGGACACCGTAGGAGACCAATACTGCGTAAATCATAAGCATTGGAGAGTGGACTGCTAGCGAGGCTCACCATGGAGTCTCGGTCGGCATCTCCTGACTGCTGTTCCATCGCGTTTTTCTTTTACTCACGCAATAAATCAATACCCCCTAACACAGGCCTGCTCCAGCCTTATTAAGGCCATAGTAGCTCTACATGTAGACCGAACGGAAGCACAGTTTGGTAGAAATTCTTAATCGACTATGGTCCGTGCAGGCCAAAAAAGGAATAATCTTCGAATTCTCACGCCTTCATTAGGGCGCACATGGTGGGGTAAATCACTGCACTCTGTTCGCAGTTAAGCGTTGCAATCAATATCGGCAGAACTCGGAGTCCGTATAAAGCCGCCTCAGCGTGCACACGCCCGTGCGGCACGTCATTAGACGAGGATTCCGGGGGACTGGCCTGTTCGTAATCCACTAAAACAATGGTCCTACCATCTAAAACGCACCGTGTTCCCCTCTACGGGAACCCCCTAGAT',
                     'AGAGCGCAGAGAAGTCATTAGAACATGTAGCACATCGCTTATTAAGGGTCAATACCTAAAGGGCCTAACTATACGCCACACGGAACAGCTC')     
      print (d)
      print (s1)
      print (t1)
# Expect:
# 22
# AGGGCGCACATG--GTGGGGTA-AATCAC-T-GCAC-TCTG-TTCGCAGTTAAGCGTTGCAATCAATATCGGC-AGAACTCGGAGTCCGTA--TAAAGCCGCCTCAGCGTGCACACGC-C
# AGAGCGCAGA-GAAGTCAT-TAGAA-CATGTAGCACATC-GCTT---A-TTAAG-G--G---TCAATA-C--CTA-AA---GG-G-CC-TAACTATA-C-GCCACA-CG-GAACA-GCTC
