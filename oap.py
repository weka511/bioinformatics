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

# OAP Overlap Alignment

from align import create_distance_matrix
from rosalind_old import highest_scoring_local_alignment
from reference_tables import createSimpleDNASubst

def oap(s,t):
    return highest_scoring_local_alignment(s,t,
                                           weights=createSimpleDNASubst(subst=2),
                                           sigma=2)
if __name__=='__main__':
    from helpers import create_strings
    d,s1,t1 = oap('CTAAGGGATTCCGGTAATTAGACAG','ATAGACCATATGTCAGTGACTGTGTAA')
    #strings  = create_strings('oap',fasta=1)
    #d,s1,t1 = overlap_assignment(strings[0],strings[1])      
    print ('{0}'.format(d))
    print (s1)
    print (t1)    
    