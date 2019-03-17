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

#    SCC Strongly Connected Component
#

def dfs(adj         = None,
        n           = None,
        sequence    = None,
        previsit    = lambda v:None,
        postvisit   = lambda v:None,
        pre_explore = lambda v:None):

  def explore(v):
    visited[v] = True
    previsit(v)

    for u in adj[v]:
      if not visited[u]:
        explore(u)
          
    postvisit(v)   
  
  visited = [False for v in range(n+1)]  # Zero element won't be used, but it does simplify indexing

  for v in sequence:
    if not visited[v]:
      pre_explore(v)
      explore(v)

    
def scc(edges):
  
  def counter():
    count = 0
    while True:
      yield count
      count+=1
    
  def create_adj(edges,reverse=False):
    n,_=edges[0]
    adj = {}
    for i in range(1,n+1):
      adj[i]=[]
    for (a,b) in edges[1:]:
      if reverse:
        (a,b)=(b,a)
      adj[a].append(b)
    return adj
    
  def decreasing(post):
    pairs = sorted(zip(post,range(1,len(post)+1)),reverse=True)
    for a,b in pairs:
      yield b

  def incr_pre(v):
    pre[v]     = next(clock)
    
  def incr_post(v):
    post[v]     = next(clock)

  def incr(v):
    ccnum[v-1]=next(cc)
    
  n,_     = edges[0]
  clock   = counter()
  pre     = [-1 for v in range(n+1)]   # Zero element won't be used, but it does simplify indexing
  post    = [-1 for v in range(n+1)]   # Zero element won't be used, but it does simplify indexing 
     
  dfs(adj      = create_adj(edges,reverse=True),
      n        = n,
      sequence = range(1,n+1),
      previsit = incr_pre,
      postvisit = incr_post)
  
  cc       = counter()
  ccnum    = [-1 for v in range(n+1)]   # Zero element won't be used, but it does simplify indexing

  dfs(adj         = create_adj(edges),
      n           = n,
      sequence    = decreasing(post[1:]),
      pre_explore = incr)
  
  return len([cc for cc in ccnum if cc>-1])     

if __name__=='__main__':
  from helpers import create_list
  #print (
    #scc([
      #(6, 7),
      #(4, 1),
      #(1, 2),
      #(2, 4),
      #(5, 6),
      #(3, 2),
      #(5, 3),
      #(3, 5)
  #]))
  
  print (scc(create_list(ext=4)))
               

