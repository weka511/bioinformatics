import MER

def ms(n,A):
    if n<2:
        return A
    else:
        n1=n//2
        A1=A[0:n1]
        B1=A[n1:]
        m1=len(B1)
        return MER.mer(n1,ms(n1,A1),m1,ms(m1,B1))
    
#print (ms(10,[20, 19, 35, -18, 17, -20, 20, 1 ,4, 4]))

with open('c:/Users/Weka/Downloads/rosalind_ms(1).txt') as f:
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
    #print (n,a)
    #print ()
    print (' '.join([str(r) for r in ms(n,a)]))