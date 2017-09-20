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

# BA8B Compute the Squared Error Distortion

import random,math

def SquaredErrorDistortion(centres,data,k,m):
   def distance(p1,p2):
      return sum([(p1[i]-p2[i])**2 for i in range(m)])
 
   
   def distance_to_nearest(pt):
      return min([distance(pt,c) for c in centres])
   
   return sum([distance_to_nearest(pt) for pt in data])/len(data)


if __name__=='__main__':
   m = -1
   k = -1
   points=[]

   with open (r'C:\Users\Weka\Downloads\rosalind_ba8b(1).txt') as f:   
      for line in f:
         if k==-1:
               values=line.strip().split()
               k=int(values[0])
               m=int(values[1])
         else:
            if line[0:2]!='--':
               points.append([float(v) for v in line.strip().split()])
   print( SquaredErrorDistortion(points[0:k],points[k:],k,m))