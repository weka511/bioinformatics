def par(n,A):
    prefix=[]
    suffix=[A[0]]
    for a in A[1:]:
        if a<=A[0]:
            prefix.append(a)
        else:
            suffix.append(a)
    return prefix + suffix

#print(par(9,[7, 2 ,5 ,6, 1, 3, 9, 4, 8]))

if __name__=='__main__':

    with open('c:/Users/Weka/Downloads/rosalind_par.txt') as f:
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
   
    print (' '.join([str(r) for r in par(n,A)]))