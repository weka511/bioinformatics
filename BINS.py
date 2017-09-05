def bins(n,m,a,values):
    def search(value,imin,imax):
        if imin>imax:
            return -1
        elif imin==imax:
            return imin if value==a[imin-1] else -1
        else:
            imid=(imax+imin)//2
            left=search(value,imin,imid)
            if left>-1:
                return left
            else:
                return search(value,imid+1,imax)
    return [search(value,1,n) for value in values]


with open('c:/Users/Weka/Downloads/rosalind_bins(1).txt') as f:
    i=0
    n=0
    m=0
    a=[]
    values=[]
    for line in f:
        text=line.strip()
        if i==0:
            n=int(text)
        elif i==1:
            m=int(text)
        elif i==2:
            a=[int(t) for t in text.split(' ')]
            print (a)
        else:
            values=[int(t) for t in text.split(' ')]
        i+=1
    indices = bins(n,m,a,values)
    
print (' '.join([str(r) for r in indices]))