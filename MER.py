def mer(n,A,m,B):
    i=0
    j=0
    C=[]
    while j<n and i<m:
        if A[j]<B[i]:
            C.append(A[j])
            j+=1
        else:
            C.append(B[i])
            i+=1
    while j<n:
        C.append(A[j])
        j+=1
    while i<m:
        C.append(B[i])
        i+=1
    return C

#print (mer(4,[2, 4, 10, 18],3,[-5 ,11, 12]))

with open('c:/Users/Weka/Downloads/rosalind_mer.txt') as f:
    n=0
    A=[]
    m=0
    B=[]
    i=0
    for line in f:
        text=line.strip()
        if i==0:
            n=int(text)
        elif i==1:
            A=[int(t) for t in text.split(' ')]
        elif i==2:
            m=int(text)
        else:
            B=[int(t) for t in text.split(' ')]
        i+=1
    
print (' '.join([str(r) for r in mer(n,A,m,B)]))