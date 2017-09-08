def med(n,A,k):
    def helper(n,A):
        a_min=A[0]
        index=0
        for kk in range(n):
            if A[kk]<a_min:
                index=kk
                a_min=A[index]
        if index==0:
            return (a_min,index,A[1:])
        if index==n-1:
            return (a_min,index,A[0:-1])
        else:
            return (a_min,index,A[0:index]+A[index+1:])
    AA=A[:]    
    for kk in range(k):
            a,_,AA=helper(len(AA),AA)
            #if kk%100==0:
                #print (kk,a)
       
    return a

#print (med(11,[2, 36, 5, 21, 8, 13, 11, 20, 5, 4 ,1],8))

if __name__=='__main__':
    import timeit
    start_time = timeit.default_timer()
    with open('c:/Users/Weka/Downloads/rosalind_med.txt') as f:
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
        print (med(n,A,k))
        print ('Elapsed: {0} seconds'.format(timeit.default_timer() - start_time))