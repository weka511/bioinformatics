# Copyright (C) 2015-2019 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# This file contains a collection of functions to solve the problems
# at rosalind.info.

import helpers as rh, reference_tables as rrt, fasta as f, math, sys, numpy,random, functools

class RosalindException(Exception):
    def __init__( self, message ):
        Exception.__init__(self, message)
        
## Rosalind functions

### How do we sequence antibodies?




	  	  	 

### Test cases ###

if __name__=='__main__':
 
    import unittest,fragile
    
    class TestRosalind(unittest.TestCase):
            




# BA5G 	Compute the Edit Distance Between Two Strings 	 	
                         
# BA5H 	Find a Highest-Scoring Fitting Alignment of Two Strings 	 	
                         
# BA5I 	Find a Highest-Scoring Overlap Alignment of Two Strings 	 	
                         
# BA5J 	Align Two Strings Using Affine Gap Penalties 	 	
                         
# BA5K 	Find a Middle Edge in an Alignment Graph in Linear Space 	 	
                         
# BA5L 	Align Two Strings Using Linear Space 	 	
                         
# BA5M 	Find a Highest-Scoring Multiple Sequence Alignment 	 	
                         
# BA5N 	Find a Topological Ordering of a DAG 
        def test_ba5n(self):
            self.assertEqual([5, 4, 1, 2, 3],
                             topological_order({
                                 1 : [2],
                                 2 : [3],
                                 4 : [2],
                                 5 : [3]
                             }))
            
#LGIS 	Longest Increasing Subsequence
        def test_longestIncreasingSubsequence(self):
            (a,d)=longestIncreasingSubsequence(5,[5, 1, 4, 2, 3])
            self.assertEqual([1,2,3],a)
            self.assertEqual([5,4,3],d)
            (a,d)=longestIncreasingSubsequence(9,[8, 2, 1, 6, 5, 7, 4, 3, 9])
            self.assertEqual([1, 5, 7, 9],a)
            self.assertEqual([8, 6, 5, 4, 3],d)
        
        def test_orf(self):
            string='''>Rosalind_99
            AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'''
            fasta=f.FastaContent(string.split('\n'))            
            peptides=get_reading_frames(fasta)
            self.assertIn('MLLGSFRLIPKETLIQVAGSSPCNLS',peptides)
            self.assertIn('M',peptides)
            self.assertIn('MGMTPRLGLESLLE',peptides)
            self.assertIn('MTPRLGLESLLE',peptides)            
            self.assertEqual(4,len(peptides))
            #for peptide in get_reading_frames(f.FastaFile('./data/rosalind_orf(3).txt') ):
                #print (peptide)

        #def test_wfmd(self):
            #self.assertAlmostEqual(0.772, wfmd(4,6,2,1),3)
            
        #def test_q5(self):
            #for peptide in ['MTAI',
                            #'MLAT',
                            #'MAIT',
                            #'MTAL',
                            #'TALM',
                            #'TLAM']:
                #sp=cycloSpectrum(peptide)
                #if sp==[0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]:
                    #print (peptide, sp)
         
        #def test_q6(self):
            ## 0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333
            #spectrum=[0 ,71, 99, 101, 103, 128, 129, 199, 200, 204, 227, 230, 231,\
                #298, 303, 328, 330, 332, 333]
            #print (spectrum)
            #for peptide in ['TVQ',
                            #'CTV',
                            #'AVQ',
                            #'AQV',
                            #'TCE',
                            #'VAQ']:
                #masses=[rrt.integer_masses[a] for a in peptide]
                #masses.sort()
                #ls=rh.linearSpectrum(masses)
                #if rh.consistent(masses,spectrum):
                    #print (peptide)
                #else:
                    #print (peptide,masses,ls)
                    
    unittest.main(exit=False)
