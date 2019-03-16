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

#    NWC Negative Weight Cycle
#

    
def scc(edges):
  
  def counter():
    count = 0
    while True:
      yield count
      count+=1
    
  def create_adj_R(edges):
    adj_R = {}
    for (a,b) in edges[1:]:
      if not b in adj_R:
        adj_R[b]=[]
      adj_R[b].append(a)
    return adj_R
 
  def dfs(adj,n,previsit=lambda v:None,postvisit=lambda v:None):
    clock = counter()
    def explore(v):
      visited[v] = True
      previsit(v)
      pre[v]     = next(clock)
      for u in adj[v+1]:
        if not visited[u]:
          explore(u-1)
      postvisit(v)   
      post[v]    = next(clock)
      return
    
    visited = [False for v in range(n)]
    pre     = [-1 for v in range(n)]
    post    = [-1 for v in range(n)]
    for v in range(n):
      if not visited[v]:
        explore(v)
    return
  
  n,_ = edges[0]
  adj_R = create_adj_R(edges)
  dfs(adj_R,n)
  return -1     

if __name__=='__main__':
  edges = [(6, 7),
           (4, 1),
           (1, 2),
           (2, 4),
           (5, 6),
           (3, 2),
           (5, 3),
           (3, 5)
           ]
  print (scc(edges))
               

