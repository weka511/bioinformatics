def step(centres,points):
    return centres

def soft_kmeans(k,m,beta,points,steps=100):
    centres=points[:k]
    for i in range(steps):
        centres=step(centres,points) 
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