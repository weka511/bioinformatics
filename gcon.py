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

# GCON Global Alignment with Constant Gap Penalty

from align import align

if __name__=='__main__':
    from Bio.SubsMat.MatrixInfo import blosum62
    score,s1,s2=align('PLEASANTLY','MEANLY',replace_score=blosum62,indel_cost=(5,0))
    print (score,s1,s2)
