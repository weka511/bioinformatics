import numpy as np

def soft_kmeans(k,m,beta,points,steps=1):
    def distance(p1,p2):
        return np.sqrt(sum((p1[i]-p2[i])**2 for i in range (m)))    
    @np.vectorize
    def hidden_matrix(i,j):
        return np.exp(-beta*distance(centres[i],points[j]))
    def step(centres):
        ks=np.array(range(k))
        ms=np.array(range(m))
        kk,mm=np.meshgrid(ks,ms)
        numerator=hidden_matrix(kk,mm)
        print (numerator)
        denominator=np.sum(numerator,axis=0)
        print (denominator)
        matrix=np.divide(numerator,denominator)
        print (matrix)
        return centres

    centres=points[:k]
    for i in range(steps):
        centres=step(centres) 
    return centres

if __name__=='__main__':
    m = -1
    k = -1
    beta = -1
    points=[]
 
    #with open (r'C:\Users\Weka\Downloads\rosalind_ba8c.txt') as f:   
    with open('BA8D.txt') as f:
        for line in f:
            if k==-1:
                values=line.strip().split()
                k=int(values[0])
                m=int(values[1])
            elif beta==-1:
                beta=float(line.strip())
            else:
                points.append([float(v) for v in line.strip().split()])
        for pt in soft_kmeans(k,m,beta,points):
            print (' '.join('{0:.3f}'.format(p) for p in pt))