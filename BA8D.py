

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