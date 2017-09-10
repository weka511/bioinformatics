def match(s1,s2):
    index=-1
    for i in range(min(len(s1),len(s2))//2,min(len(s1),len(s2))):
        if s1[len(s1)-i:]==s2[:i]:
            index=i
    return index

def long(strings):
    def compare():
        pairs=[]
        for i in range(len(strings)):
            for j in range(len(strings)):
                if i!=j:
                    _,s1=strings[i]
                    _,s2=strings[j]
                    m=match(s1,s2)
                    if m>0:
                        pairs.append((i,j,m,s1,s2))
        return pairs
    def split(pairs):
        froms=[]
        tos=[]
        for (a,b,_,_,_) in pairs:
            froms.append(a)
            tos.append(b)
        return (froms,tos)
    def get_unmatched(froms,tos):
        for i in range(len(froms)):
            matched=False
            for j in range(len(tos)):
                if tos[j]==froms[i]:
                    matched=True
                    break
            if not matched:
                return i
            
    pairs=compare()
    genome=[]
    while len(pairs)>0:
        (froms,tos)=split(pairs)
        index=get_unmatched(froms,tos)
        pair=pairs[index]
        _,_,length,pre,post=pair
        del pairs[index]
        if len(genome)==0:
            genome.append(pre)
        genome.append(post[length:])
    return ''.join(genome)

if __name__=='__main__':
    def combine(seq_record):
        return (seq_record.id,str(seq_record.seq))
    
    from Bio import SeqIO
    print (long([combine(seq_record) for seq_record in SeqIO.parse("LONG.txt", "fasta")]))
