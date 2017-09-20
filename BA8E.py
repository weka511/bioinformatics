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

# BA8E Implement Hierarchical Clustering

def HierarchicalClustering(D, n):
  pass


if __name__=='__main__':
  print (HierarchicalClustering(
   [[0.00,0.74,0.85,0.54,0.83,0.92,0.89],
    [0.74,0.00,1.59,1.35,1.20,1.48,1.55],
    [0.85,1.59,0.00,0.63,1.13,0.69,0.73],
    [0.54,1.35,0.63,0.00,0.66,0.43,0.88],
    [0.83,1.20,1.13,0.66,0.00,0.72,0.55],
    [0.92,1.48,0.69,0.43,0.72,0.00,0.80],
    [0.89,1.55,0.73,0.88,0.55,0.80,0.00]],
   7))
 