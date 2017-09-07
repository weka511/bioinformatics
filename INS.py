def ins(n,A):
    swap=0
    for i in range(1,n):
        k=i
        while k>0 and A[k]<A[k-1]:
            A[k-1],A[k]=A[k],A[k-1]
            k-=1
            swap+=1

    return (swap,A)

with open('c:/Users/Weka/Downloads/rosalind_ins.txt') as f:
    i=0
    n=0
    a=[]
    for line in f:
        text=line.strip()
        if i==0:
            n=int(text)
        else:
            a=[int(t) for t in text.split(' ')]
        i+=1
    print (ins(n,a))
