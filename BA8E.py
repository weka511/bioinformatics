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

def minimum_distance(Cluster1,Cluster2):
  if len(Cluster1)==0 or len(Cluster2)==0:
    return float('inf')
  return min([distance(Cluster1[i],Cluster2[j]) for i in range(len(Cluster1)) for j in range(len(Cluster2))  ])

def average_distance(Cluster1,Cluster2):
  if len(Cluster1)==0 or len(Cluster2)==0:
    return float('inf')
  return sum([distance(Cluster1[i],Cluster2[j]) for i in range(len(Cluster1)) for j in range(len(Cluster2))  ])/(len(Cluster1)*len(Cluster2))
  
def HierarchicalClustering(D, n,cluster_distance=average_distance):
  def create_graph():
    T = {}
    for i in range(n):
      T[i]=[]
    return T
  

  
  def find_two_closest_clusters(Clusters,cluster_distance=minimum_distance):
    m=-1
    n=-1
    distance=float('inf')
    for i in range(len(Clusters)):
      for j in range(i):
        d=cluster_distance(Clusters[i],Clusters[j])
        if d<distance:
          m=i
          n=j
          distance=d
    return (m,n,distance)
  
  Clusters = [[point] for point in D]
  T= create_graph()
  
  for k in range(n-1):
    i,j,distance=find_two_closest_clusters(Clusters) #i>j
    C_new=Clusters[i]+Clusters[j]
    Clusters.append(C_new)
    T[len(T)]=[i,j]
    Clusters[i]=[]
    Clusters[j]=[]  
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
   7,cluster_distance=minimum_distance):
    print (e)
    
    
    #for e in HierarchicalClustering(
      #[[0.00,0.43,0.93,0.85,0.94,0.70,0.95,0.79,0.89,1.34,0.80,0.64,1.42,1.37,0.86,1.20,0.49,0.19,1.48,0.86],
      #[0.43,0.00,0.66,0.70,1.48,0.58,0.70,0.55,0.68,1.45,1.31,0.37,1.24,1.69,0.91,1.28,1.04,0.29,0.94,1.05],
      #[0.93,0.66,0.00,1.01,1.42,0.80,0.76,1.18,0.21,1.18,1.35,0.82,0.68,0.85,0.75,1.01,0.93,1.03,0.60,0.91],
      #[0.85,0.70,1.01,0.00,1.02,0.22,0.89,0.91,0.79,1.47,1.04,1.51,0.71,1.63,0.34,1.05,1.24,0.56,1.09,1.52],
      #[0.94,1.48,1.42,1.02,0.00,1.46,1.48,0.95,1.44,1.08,0.30,1.29,1.28,0.65,1.00,0.80,0.80,0.92,1.48,0.65],
      #[0.70,0.58,0.80,0.22,1.46,0.00,0.77,1.00,0.65,1.25,1.27,1.36,0.81,1.56,0.37,0.98,1.17,0.68,0.88,1.23],
      #[0.95,0.70,0.76,0.89,1.48,0.77,0.00,1.41,1.08,1.61,1.65,0.89,0.76,1.32,0.64,0.70,0.59,1.07,0.49,0.91],
      #[0.79,0.55,1.18,0.91,0.95,1.00,1.41,0.00,1.05,0.73,1.08,0.70,1.49,1.17,1.00,1.34,1.45,0.49,1.06,1.30],
      #[0.89,0.68,0.21,0.79,1.44,0.65,1.08,1.05,0.00,0.96,1.09,0.94,0.44,1.06,0.90,1.47,1.20,0.79,1.04,1.39],
      #[1.34,1.45,1.18,1.47,1.08,1.25,1.61,0.73,0.96,0.00,0.96,1.10,1.05,0.48,1.36,1.26,1.38,1.38,1.03,1.38],
      #[0.80,1.31,1.35,1.04,0.30,1.27,1.65,1.08,1.09,0.96,0.00,1.08,1.09,0.79,1.40,1.03,1.02,0.78,1.79,0.86],
      #[0.64,0.37,0.82,1.51,1.29,1.36,0.89,0.70,0.94,1.10,1.08,0.00,1.42,1.20,1.61,1.30,0.86,0.68,1.04,0.83],
      #[1.42,1.24,0.68,0.71,1.28,0.81,0.76,1.49,0.44,1.05,1.09,1.42,0.00,0.99,0.84,1.20,1.21,1.22,0.97,1.58],
      #[1.37,1.69,0.85,1.63,0.65,1.56,1.32,1.17,1.06,0.48,0.79,1.20,0.99,0.00,1.13,0.61,1.00,1.60,0.81,0.83],
      #[0.86,0.91,0.75,0.34,1.00,0.37,0.64,1.00,0.90,1.36,1.40,1.61,0.84,1.13,0.00,0.68,0.89,0.86,0.65,1.04],
      #[1.20,1.28,1.01,1.05,0.80,0.98,0.70,1.34,1.47,1.26,1.03,1.30,1.20,0.61,0.68,0.00,0.92,1.47,0.47,0.42],
      #[0.49,1.04,0.93,1.24,0.80,1.17,0.59,1.45,1.20,1.38,1.02,0.86,1.21,1.00,0.89,0.92,0.00,0.90,1.18,0.49],
      #[0.19,0.29,1.03,0.56,0.92,0.68,1.07,0.49,0.79,1.38,0.78,0.68,1.22,1.60,0.86,1.47,0.90,0.00,1.56,1.21],
      #[1.48,0.94,0.60,1.09,1.48,0.88,0.49,1.06,1.04,1.03,1.79,1.04,0.97,0.81,0.65,0.47,1.18,1.56,0.00,0.84],
      #[0.86,1.05,0.91,1.25,0.65,1.23,0.91,1.30,1.39,1.38,0.86,0.83,1.58,0.83,1.04,0.42,0.49,1.21,0.84,0.00] ]     ,
     #20,cluster_distance=average_distance):
      #print (e)
  
  #with open (r'C:\Users\Weka\Downloads\rosalind_ba8e.txt') as f: 
    #n=-1
    #D=[]
    #for line in f:
      #if n==-1:
        #n=int(line.strip())
      #else:
          #D.append([float(v) for v in line.strip().split()])  
  #for e in HierarchicalClustering(D, n,cluster_distance=average_distance):
    #print (e)
 