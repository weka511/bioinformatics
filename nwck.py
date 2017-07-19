import newick

tokenizer = newick.Tokenizer()
parser = newick.Parser(tokenizer)

def get_ancestors(clade,path=[]):
    ancestor = clade.parent
    if ancestor == None:
        return path
    pp=path[:]
    pp.append(ancestor.id)
    return get_ancestors(ancestor,pp)

def get_path(clade):
    return [clade.id]+get_ancestors(clade)

def diff(path1,path2):
    if len(path1)>len(path2):
        return diff(path2,path1)
    
    i=0
    while i<len(path1)  and path1[i]==path2[i]:
        i+=1
 
    return len(path1[i:]) +len(path2[i:])
    
with open (r'C:\Users\Weka\Downloads\rosalind_nwck(1).txt') as f:
    diffs=[]
    i = 0
    tree = None
    lookup = None
    for line in f:
        if i%3==0:
            print (line.strip())
            tree,lookup=parser.parse(line.strip())
        elif i%3==1:
            clades = line.strip().split()
            print (clades)
            ancestors=[get_path(lookup[clade])[::-1] for clade in clades]
            diffs.append(diff(ancestors[0],ancestors[1]))
        i+=1
    print (diffs)