
import random

def selection(A,k):
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
    v=random.choice(A)
    (prefix,mid,suffix)=split(v)
    if k<=len(prefix):
        return selection(prefix,k)
    elif len(prefix)<k and k<=len(prefix)+len(mid):
        return v
    else:
        return selection(suffix, k-len(prefix)-len(mid))
    
#print (selection([2, 36, 5, 21, 8, 13, 11, 20, 5, 4 ,1],8))

if __name__=='__main__':
    import timeit
    start_time = timeit.default_timer()
    with open('c:/Users/Weka/Downloads/rosalind_med(1).txt') as f:
        n=0
        A=[]
        k=0
        i=0
        for line in f:
            text=line.strip()
            if i==0:
                n=int(text)
            elif i==1:
                A=[int(t) for t in text.split(' ')]
            elif i==2:
                k=int(text)
            i+=1
        print('n={0},k={1}'.format(n,k))
        print (selection(A,k))
        print ('Elapsed: {0} seconds'.format(timeit.default_timer() - start_time))