import numpy as np

def soft_kmeans(k,m,beta,points,N=1):
    def distance(p1,p2):
        return np.sqrt(sum((p1[i]-p2[i])**2 for i in range (m)))    
    @np.vectorize
    def hidden_matrix(i,j):
        return np.exp(-beta*distance(centres[i],points[j]))
    def step(centres):
        ks=np.array(range(k))
        ns=np.array(range((len(points))))
        kk,mm=np.meshgrid(ns,ks,indexing='xy')
        numerator=hidden_matrix(mm,kk)
        denominator=np.sum(numerator,axis=0)
        matrix=np.divide(numerator,denominator)
        #print (matrix)
        new_centres=[[] for i in range(k)]
        for i in range(k):
            for j in range(m):
                x_i_j=sum(matrix[i,l]*points[l][j] for l in range(len(points)))/sum(matrix[i,l] for l in range(len(points)))
                new_centres[i].append(x_i_j)
        return new_centres

    centres=points[:k]
    for i in range(N):
        centres=step(centres) 
    return centres

if __name__=='__main__':
    m = -1
    k = -1
    beta = -1
    points=[]
 
    with open (r'C:\Users\Weka\Downloads\rosalind_ba8d.txt') as f:   
    #with open('BA8D.txt') as f:
        for line in f:
            if k==-1:
                values=line.strip().split()
                k=int(values[0])
                m=int(values[1])
            elif beta==-1:
                beta=float(line.strip())
            else:
                points.append([float(v) for v in line.strip().split()])
        for pt in soft_kmeans(k,m,beta,points,N=100):
            print (' '.join('{0:.3f}'.format(p) for p in pt))