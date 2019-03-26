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

from align import FindHighestScoringMultipleSequenceAlignment




if __name__=='__main__':
    from helpers import create_strings    
    s,u,v,w = FindHighestScoringMultipleSequenceAlignment('ATATCCG','TCCGA','ATGTACTG')
    print (s)
    print (u)
    print (v)
    print (w)
    
    s1,u1,v1,w1 = FindHighestScoringMultipleSequenceAlignment('TGTTTAAAAATGTCCGCAACCATTTC',
                                                              'GATATAAAACAGGGATAACTGCAATGG',
                                                              'CCTGCTACTTTATGCCGTCTCCATATGCG')
    print (s1)
    print (u1)
    print (v1)
    print (w1) 
    
    ss=create_strings(ext=3)
    s2,u2,v2,w2 = FindHighestScoringMultipleSequenceAlignment(ss[0],ss[1],ss[2])
    print (s2)
    print (u2)
    print (v2)
    print (w2)    
 