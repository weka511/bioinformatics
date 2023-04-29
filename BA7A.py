#!/usr/bin/env python

# Copyright (C) 2017-2023 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

'''BA7A Compute Distances Between Leaves'''

from helpers   import create_weighted_adjacency_list
from phylogeny import ComputeDistancesBetweenLeaves


if __name__=='__main__':

    n,T = create_weighted_adjacency_list()
    for ds in ComputeDistancesBetweenLeaves(n,T):
        print(' '.join(str(d) for d in ds))
