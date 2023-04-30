#!/usr/bin/env python

# Copyright (C) 2017-2023 Greenweaves Software Pty Ltd

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



from phylogeny import NeighborJoining
from rosalind  import Tree,read_matrix


if __name__=='__main__':
    params,D=read_matrix('c:/Users/Weka/Downloads/rosalind_ba7e.txt',conv=float)
    NeighborJoining(D,params[0]).print()
