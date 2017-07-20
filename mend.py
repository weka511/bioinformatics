import newick

def mend(node):
    if len(node.nodes)==0:
        print (node.name)
    else:
        for parent in node.nodes:
            mend(parent)
            
tokenizer = newick.Tokenizer()
parser = newick.Parser(tokenizer)

tree,_=parser.parse('((((Aa,aa),(Aa,Aa)),((aa,aa),(aa,AA))),Aa);')
mend(tree)

