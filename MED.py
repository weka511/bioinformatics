def med(n,A,k):
    #print (n,A,k)
    if k==1:
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
        
    else:
        (a_min,index,A_rest)=med(n-1,A,k-1)
        return med(len(A_rest),A_rest,1)
    
#print (med(11,[2, 36, 5, 21, 8, 13, 11, 20, 5, 4 ,1],8))

if __name__=='__main__':

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
        print (med(n,A,k))