import random

def qs(A):
    def split(v):
        prefix=[]
        suffix=[]
        mid=[]
        for a in A:
            if a<v:
                prefix.append(a)
            elif a==v:
                mid.append(a)
            else:
                suffix.append(a)
        return (prefix,mid,suffix)

    if len(A)<2:
        return A
    v=random.choice(A)
    (prefix,mid,suffix)=split(v)
    return qs(prefix)+mid+qs(suffix)

#print (qs([5, -2, 4, 7, 8, -10, 11]))

if __name__=='__main__':
    import timeit
    start_time = timeit.default_timer()
    with open('c:/Users/Weka/Downloads/rosalind_qs.txt') as f:
        n=0
        A=[]
        i=0
        for line in f:
            text=line.strip()
            if i==0:
                n=int(text)
            elif i==1:
                A=[int(t) for t in text.split(' ')]
            i+=1
        print (qs(A))
        print ('Elapsed: {0} seconds'.format(timeit.default_timer() - start_time))