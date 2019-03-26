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

# BA5M.py Find a Highest-Scoring Multiple Sequence Alignment 


from Bio.SubsMat.MatrixInfo import blosum62


def FindHighestScoringMultipleSequenceAlignment (v,w,x,replace_score=blosum62,indel_cost=5):
    pass

if __name__=='__main__':
    from helpers import create_strings    
    print (FindHighestScoringMultipleSequenceAlignment('ATATCCG','TCCGA','ATGTACTG'))
 