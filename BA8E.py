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

from rosalind import distance

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

def cdistance(Cluster1,Cluster2):
  if len(Cluster1)==0 or len(Cluster2)==0:
    return float('inf')
  return min([distance(Cluster1[i],Cluster2[j]) for i in range(len(Cluster1)) for j in range(len(Cluster2))  ])

def udistance(Cluster1,Cluster2):
  if len(Cluster1)==0 or len(Cluster2)==0:
    return float('inf')
  return sum([distance(Cluster1[i],Cluster2[j]) for i in range(len(Cluster1)) for j in range(len(Cluster2))  ])/(len(Cluster1)*len(Cluster2))
  
def HierarchicalClustering(D, n,cluster_distance=udistance):
  def create_graph():
    T = {}
    for i in range(n):
      T[i]=[]
    return T
  

  
  def find_two_closest_clusters(Clusters,cluster_distance=cdistance):
    m=-1
    n=-1
    minimum_distance=float('inf')
    for i in range(len(Clusters)):
      for j in range(i):
        d=cluster_distance(Clusters[i],Clusters[j])
        if d<minimum_distance:
          m=i
          n=j
          minimum_distance=d
    return (m,n,minimum_distance)
  
  Clusters = [[point] for point in D]
  T= create_graph()
  
  for k in range(n-1):
    i,j,minimum_distance=find_two_closest_clusters(Clusters) #i>j
    C_new=Clusters[i]+Clusters[j]
    Clusters.append(C_new)
    T[len(T)]=[i,j]
    Clusters[i]=[]
    Clusters[j]=[]
    for c in Clusters:
      print (c)
  print ("Done")
  return (format(T))
  
if __name__=='__main__':
  for e in HierarchicalClustering(
   [[0.00,0.74,0.85,0.54,0.83,0.92,0.89],
    [0.74,0.00,1.59,1.35,1.20,1.48,1.55],
    [0.85,1.59,0.00,0.63,1.13,0.69,0.73],
    [0.54,1.35,0.63,0.00,0.66,0.43,0.88],
    [0.83,1.20,1.13,0.66,0.00,0.72,0.55],
    [0.92,1.48,0.69,0.43,0.72,0.00,0.80],
    [0.89,1.55,0.73,0.88,0.55,0.80,0.00]],
   7):
    print (e)
  
 