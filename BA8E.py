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

from rosalind import distance,read_matrix,write_list

def format(T):
  def expand(links):
    result=[]
    for link in links:
      ll=T[link]
      if len(ll)==0:
        result.append(link+1)
      else:
        result=result+expand(ll)
    result.sort()
    return result
  
  return [expand(links) for key,links in T.items() if len(links)>0]

  
def HierarchicalClustering(D, n):
  def create_graph():
    T = {}
    for i in range(n):
      T[i]=[]
    return T
 
  def find_two_closest_clusters():
    ii=-1
    jj=-1
    best_distance=float('inf')
    for i in range(len(D)):
      for j in range(i):
        if D[i][j]<best_distance:
          ii=i
          jj=j
          best_distance=D[i][j]
    return (ii,jj,best_distance)
  
  NClusters=len(D)
  T= create_graph()
  
  for k in range(n-1):
    def d(i,j,k):
      return (D[i][k]+D[j][k])/2
    i,j,distance=find_two_closest_clusters() #i>j
    NClusters+=1
    T[len(T)]=[i,j]

    row=[d(i,j,k) for k in range(len(D))] + [0.0]
    for k in range(len(D)):
      D[k].append(row[k])
    D.append(row)
    D[i]=[float('inf')]*len(D)
    D[j]=[float('inf')]*len(D) 
    for ds in D:
      ds[i]=float('inf')
      ds[j]=float('inf')
  return (format(T))
  
if __name__=='__main__':
 
  params,D=read_matrix('c:/Users/Weka/Downloads/rosalind_ba8e(4).txt',conv=float) 

  for e in HierarchicalClustering(D, params[0]):
    write_list(e)
 