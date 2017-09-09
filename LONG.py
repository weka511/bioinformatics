def match(s1,s2):
    index=-1
    for i in range(min(len(s1),len(s2))//2,min(len(s1),len(s2))):
        if s1[len(s1)-i:]==s2[:i]:
            index=i
    #if index>-1:
        #print(s1[len(s1)-index:],s2[:index])
    return index

def long(strings):
    candidate=[]
    for i in range(len(strings)):
        for j in range(len(strings)):
            if i!=j:
                _,s1=strings[i]
                _,s2=strings[j]
                m=match(s1,s2)
                if m>0:
                    candidate.append((i,j,m,s1,s2))
    print (candidate)
     
if __name__=='__main__':
    def combine(seq_record):
        return (seq_record.id,str(seq_record.seq))
    
    from Bio import SeqIO
    print (long([combine(seq_record) for seq_record in SeqIO.parse("LONG.txt", "fasta")]))
