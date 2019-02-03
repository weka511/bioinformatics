import sys

def align(x,
          y,
          s=[[ 2, -7, -5, -7],
             [-7,  2, -7, -5],
             [-5, -7,  2, -7],
             [-7, -5, -7,  2]],
          d=-5,
          bound=-sys.maxsize):
    def traceback(F):
        i    = len(x)
        j    = len(y)
        path = []
        while i > 0 or j > 0:
            predecessors = []
            indices      = []
            xs           = []
            ys           = []
            if i>0:
                predecessors.append(F[i-1][j]+d)
                indices.append((i-1,j))
                xs.append(x[i-1])
                ys.append('-')
            if j>0:
                predecessors.append(F[i][j-1]+d)
                indices.append((i,j-1))
                xs.append('-')
                ys.append(y[j-1])                
            if i>0 and j>0:
                xi          = 'ACGT'.index(x[i-1])
                yj          = 'ACGT'.index(y[j-1])                
                predecessors.append(F[i-i][j-1]+s[xi][yj])
                indices.append((i-1,j-1))
                xs.append(x[i-1])
                ys.append(y[j-1])                
            m   = predecessors.index(max(predecessors))
            i,j = indices[m]
            path.append((i,j,xs[m],ys[m]))
        return path
    
    F = [[0 for j in range(len(y)+1)] for i in range(len(x)+1)]
    for i in range(1,len(x)):
        F[i][0] = F[i-1][0] + d
    for j in range(1,len(y)):
        F[0][j] = F[0][j-1] + d    
    for i in range(len(x)):
        for j in range(len(y)):
            xi          = 'ACGT'.index(x[i])
            yj          = 'ACGT'.index(y[j])
            F[i+1][j+1] = max(F[i][j] + s[xi][yj],
                              F[i][j+1] + d,
                              F[i+1][j] + d,
                              bound)
    for f in F:
        print (f)
    t = traceback(F)
    for (_,_,x,y) in t:
        print (x,y)
    return (F[len(x)+1][len(y)+1],''.join([x for (_,_,x,_) in t]),''.join([y for (_,_,Y,y) in t]))

if __name__=='__main__':
    print (align('AGCT','AAGT'))
