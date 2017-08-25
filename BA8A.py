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
    centres = [random.choice(data)]
    while len(centres)<k:
        centres.append(furthest_point())
    return centres

if __name__=='__main__':
    print (FarthestFirstTraversal([(0.0, 0.0),
                                   (5.0, 5.0),
                                   (0.0, 5.0),
                                   (1.0, 1.0),
                                   (2.0, 2.0),
                                   (3.0, 3.0),
                                   (1.0, 2.0)
                                   ],3,2))