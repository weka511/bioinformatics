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

# BA7C Implement Additive Phylogeny

def AdditivePhylogeny(D,n):
    pass

if __name__=='__main__':
    n=4
    D=[[0,   13,  21,  22],
        [13,  0 ,  12,  13],
        [21,  12 , 0 ,  13],
        [22 , 13 , 13,  0]]
    print (AdditivePhylogeny(D,n))
