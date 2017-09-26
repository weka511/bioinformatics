# Copyright (C) 2017 Greenweaves Software Pty Ltd

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

# BA7B Limb Length Problem

from rosalind import read_matrix

def ComputeLimbLength(n,j,D):
    return int(min([D[i][j]+D[j][k]-D[i][k] for i in range(n) for k in range(n) if j!=k and k!=i and i!=j])/2)

if __name__=='__main__':
    params,D=read_matrix('c:/Users/Weka/Downloads/rosalind_ba7b.txt',len_params=2)  
    print (ComputeLimbLength(params[0],params[1],D))
