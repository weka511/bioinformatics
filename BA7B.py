# Copyright (C) 2017-2019 Greenweaves Software Limited

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

# ComputeLimbLength
#
# Inputs: n An integer n
#         j an integer between 0 and n - 1,
#         D a space-separated additive distance matrix D (whose elements are integers).
#
# Return: The limb length of the leaf in Tree(D) corresponding to row j of this distance matrix (use 0-based indexing).
#
# Uses the Limb Length Theorem: LimbLength(j) = min(D[i][j] + D[j][k]-D[i][k])/2 over all leaves i and k

def ComputeLimbLength(n,j,D):
    return int(min([D[i][j]+D[j][k]-D[i][k] for i in range(n) for k in range(n) if j!=k and k!=i and i!=j])/2)

if __name__=='__main__':
    from helpers import read_matrix
    params,D=read_matrix(len_params=2)  
    print (ComputeLimbLength(params[0],params[1],D))
