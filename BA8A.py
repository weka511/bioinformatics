import random,math

def FarthestFirstTraversal(data,k,m):
    def d(pt1,pt2):
        return math.sqrt(sum((p-q)*(p-q) for (p,q) in zip(pt1,pt2)))
    def dist(point):
        return min([d(point,c) for c in centres ])
    def furthest_point():
        best_distance = -1
        best_point = None
        for point in data:
            if dist(point)>best_distance:
                best_distance = dist(point)
                best_point = point
        return best_point
    centres = [data[0]]
   
    while len(centres)<k:
        centres.append(furthest_point())
    return centres

if __name__=='__main__':
    m = -1
    k = -1
    points=[]
 
    with open (r'C:\Users\Weka\Downloads\rosalind_ba8a(1).txt') as f:   
#    with open (r'ba8a.txt') as f:   
        for line in f:
            if k==-1:
                values=line.strip().split()
                k=int(values[0])
                m=int(values[1])
            else:
                points.append([float(v) for v in line.strip().split()])

    for pt in FarthestFirstTraversal(points,k,m):
        print (' '.join('{0:.3f}'.format(p) for p in pt))