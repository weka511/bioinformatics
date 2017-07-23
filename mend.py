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
    #print("combine ",f1,f2)
    def total(i,contribs):
        return sum(c[i] for c in contribs)
    contribs=[]
    for i in range(3):
        for j in range(3):
            #print ("F",factors[i][j],f1[i]*f2[j])
            contribs.append([f*f1[i]*f2[j] for f in factors[i][j]])
    return [total(k,contribs) for k in range(3)]

def mend(node):
    if len(node.nodes)==0:
        try:
            if node.name=='aA':
                node.name=node.name[::-1]
            freqs=frequencies[node.name]
            #print ("NN",node.name,freqs)
            return freqs
        except KeyError:
            return (0,0)
    parent_freqs = [mend(parent) for parent in node.nodes]
    #print("PF",parent_freqs)
    parent_freqs=[pp for pp in parent_freqs if len(pp)==3]
    combined= combine(parent_freqs[0],parent_freqs[1])
    #print ("C",combined)
    return combined

tokenizer = newick.Tokenizer()
parser = newick.Parser(tokenizer)

tree,_=parser.parse('((((Aa,aa),(Aa,Aa)),((aa,aa),(aa,AA))),Aa);')
print (mend(tree))

#def mend(node):
    #if len(node.nodes)==0:
        #try:
            #return frequencies[node.name]
        #except KeyError:
            #return (0,0)
    #freqs = [mend(parent) for parent in node.nodes]
    #print (freqs)
    #result= (freqs[0][0]+freqs[1][0],freqs[0][1]+freqs[1][1])    
    #return (result[0]/sum(result),result[1]/sum(result) )

