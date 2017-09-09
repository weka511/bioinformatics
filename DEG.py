def deg(n,m,A):
    degrees=[0 for a in range(n)]
    for (a,b) in A:
        print (a,b)
        degrees[a-1]+=1
        degrees[b-1]+=1
    return [d for d in degrees if d>0]

if __name__=='__main__':
    import timeit
    start_time = timeit.default_timer()
    with open('c:/Users/Weka/Downloads/rosalind_deg(1).txt') as f:
        A=[]
        for line in f:
            text=line.strip()
            pair=text.split(' ')
            print (pair)
            A.append((int(pair[0]),int(pair[1])))
        (n,m)=A[0]
        #print(A[1:])
        print (deg(n,m,A[1:]))
        print ('Elapsed: {0} seconds'.format(timeit.default_timer() - start_time))
#A=[
   #[1, 2],
   #[2, 3],
   #[6, 3],
   #[5, 6],
   #[2, 5],
   #[2, 4],
   #[4, 1]] 

#print (deg(6,7,A))

#6 7
#1 2
#2 3
#6 3
#5 6
#2 5
#2 4
#4 1

#2 4 2 2 2 2
