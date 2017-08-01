def trivial_centres(k,m,data):
    return [data[i] for i in range(k)]

def distance(p1,p2,m):
    return sum((p1[i]-p2[i])**2 for i in range (m))

def kmeans(k,m,data,centres=trivial_centres):
    centroids=centres(k,m,data)
    step_size=float("inf")
    while step_size>0:
        new_centroids=step(k,m,data,centroids)
        step_size=sum([distance(new_centroids[i],centroids[i],m) for i in range(k)] )
        centroids=new_centroids[:]
    return centroids

def step(k,m,data,centres):
    def nearest(p):
        index=-1
        best=float('inf')
        for i in range(k):
            d1=distance(p,centres[i],m)
            if d1<best:
                index=i
                best=d1
        return index
    
    def centroid(i,indices):
        #print ('centroid', i, indices)
        count=0
        point=[0 for j in range(m)]
        for j in range(len(data)):
            if indices[j]==i:
                count+=1
                for l in range(m):
                    point[l]+=data[j][l]
        return [p/max(count,1) for p in point]
    
    indices = [nearest(p) for p in data]
    
    return [centroid(i,indices) for i in range(k)]

if __name__=='__main__':
    m = -1
    k = -1
    points=[]
 
    with open (r'C:\Users\Weka\Downloads\rosalind_ba8c.txt') as f:   
    #with open('BA8C.txt') as f:
        for line in f:
            if k==-1:
                values=line.strip().split()
                k=int(values[0])
                m=int(values[1])
            else:
                points.append([float(v) for v in line.strip().split()])
        for pt in kmeans(k,m,points):
            print (' '.join('{0:.3f}'.format(p) for p in pt))