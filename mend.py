import newick

frequencies={
    'aa':(0,0,1),
    'Aa': (0,1,0),
    'AA': (1,0,0)
}

factors=[
   [[1,0,0], [0.5,0.5,0], [0,1,0]],
   [[0.5,0.5,0],[0.25, 0.5, 0.25],[0, 0.5, 0.5]],
   [[0, 1,0],[0,0.5,0.5],[0,0,1]]
]

def combine(f1,f2):
    def total(i,contribs):
        return sum(c[i] for c in contribs)
    contribs=[]
    for i in range(3):
        for j in range(3):
            contribs.append([f*f1[i]*f2[j] for f in factors[i][j]])
    return [total(k,contribs) for k in range(3)]

def mend(node):
    if len(node.nodes)==0:
        try:
            if node.name=='aA':
                node.name=node.name[::-1]
            freqs=frequencies[node.name]
            return freqs
        except KeyError:
            return (0,0)
    parent_freqs = [mend(parent) for parent in node.nodes]
    parent_freqs=[pp for pp in parent_freqs if len(pp)==3]
    return combine(parent_freqs[0],parent_freqs[1])

tokenizer = newick.Tokenizer()
parser = newick.Parser(tokenizer)

with open (r'C:\Users\Weka\Downloads\rosalind_mend.txt') as f:
    for line in f:
        tree,_=parser.parse(line.strip())
        print (mend(tree))



