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

#    BF Bellman-Ford Algorithm
#

def bf(edges,s=1):
  
  n,_         = edges[0]
  dist        = [float('inf') for i in range(1,n+1)]
  previous    = [None for i in range(1,n+1)]
  
  dist[s-1]=0
  for i in range(n-1):
    for u,v,w in edges[1:]:
      dist[v-1] = min(dist[v-1],dist[u-1]+w)
    
  return [d if d<float('inf') else 'x' for d in dist]
        

if __name__=='__main__':
  from helpers import create_list    

  print(' '.join([str(i) for i in bf(create_list())]))
  
  #edges = [(9, 13),
           #(1, 2, 10),
           #(3, 2, 1),
           #(3, 4, 1),
           #(4, 5, 3),
           #(5, 6, -1),
           #(7, 6, -1),
           #(8, 7, 1),
           #(1, 8, 8),
           #(7, 2, -4),
           #(2, 6, 2),
           #(6, 3, -2),
           #(9, 5 ,-10),
           #(9, 4, 7)]
  
  #print (' '.join(str(x) for x in bf(edges) ))
