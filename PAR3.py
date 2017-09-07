def par3(n,A):
    prefix=[]
    suffix=[]
    mid=[]
    for a in A:
        if a<A[0]:
            prefix.append(a)
        elif a==A[0]:
            mid.append(a)
        else:
            suffix.append(a)
    return prefix + mid + suffix

#print(par3(9,[4 ,5, 6, 4 ,1, 2, 5, 7, 4]))

if __name__=='__main__':

    with open('c:/Users/Weka/Downloads/rosalind_par3.txt') as f:
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
   
    print (' '.join([str(r) for r in par3(n,A)]))