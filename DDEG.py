import DEG

def ddeg(n,M,A):
    lookup=DEG.deg(n,m,A)
    sums=[0 for a in range(n)]
    for (a,b) in A:
        sums[a-1]+=lookup[b-1]
        sums[b-1]+=lookup[a-1]    
    return sums

#n=5
#m=4
#A=[[1, 2],
   #[2, 3],
   #[4, 3],
   #[2, 4]]

#print(ddeg(n,m,A))

if __name__=='__main__':
    import timeit
    start_time = timeit.default_timer()
    with open('c:/Users/Weka/Downloads/rosalind_ddeg.txt') as f:
        A=[]
        for line in f:
            text=line.strip()
            pair=text.split(' ')
            print (pair)
            A.append((int(pair[0]),int(pair[1])))
        (n,m)=A[0]
        #print(A[1:])
        print (ddeg(n,m,A[1:]))
        print ('Elapsed: {0} seconds'.format(timeit.default_timer() - start_time))