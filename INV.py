import timeit

def inv(n,A):
    count=0
    for i in range(n):
        for j in range(i+1,n):
            if A[i]>A[j]:
                count+=1
    return count

#print (inv(5,[-6, 1, 15, 8, 10]))

with open('c:/Users/Weka/Downloads/rosalind_inv.txt') as f:
    start_time = timeit.default_timer()
    i=0
    n=0
    A=[]
    for line in f:
        text=line.strip()
        if i==0:
            n=int(text)
        else:
            A=[int(t) for t in text.split(' ')]
        i+=1
    nn=inv(n,A)
    elapsed = timeit.default_timer() - start_time
    print (nn)
    print (elapsed)