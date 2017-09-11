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

import random,math

def SquaredErrorDistortion(data,k,m):
   pass

if __name__=='__main__':
   m = -1
   k = -1
   points=[]

#   with open (r'C:\Users\Weka\Downloads\rosalind_ba8a(1).txt') as f:   
   with open (r'ba8a.txt') as f:   
      for line in f:
         if k==-1:
               values=line.strip().split()
               k=int(values[0])
               m=int(values[1])
         else:
               points.append([float(v) for v in line.strip().split()])

   for pt in FarthestFirstTraversal(points,k,m):
      print (' '.join('{0:.3f}'.format(p) for p in pt))