#   Grammar from https://en.wikipedia.org/wiki/Newick_format
#
#   Tree -> Subtree ";" | Branch ";"
#   Subtree -> Leaf | Internal
#   Leaf -> Name
#   Internal -> "(" BranchSet ")" Name
#   BranchSet-> Branch | Branch "," BranchSet
#   Branch -> Subtree Length
#   Name -> empty | string
#   Length -> empty | ":" number

def parse(string):
    '''
    Tree -> Subtree ";" | Branch ";
    '''
    def parse_subtree(start,end):
        '''
         Subtree -> Leaf | Internal
        '''
        subtree,rest=parse_leaf(start,end)
        if subtree==None:
            subtree,rest=parse_internal(start,end)
        if subtree==None:
            return (None,start)
        return (subtree,rest)
    
    def parse_leaf(start,end):
        '''
         Leaf -> Name
         '''
        leaf,rest=parse_name(start,end)
        if leaf==None:
            return (None,start)        
        return (leaf,rest)
    
    def parse_internal(start,end):
        '''
        Internal -> "(" BranchSet ")" Name
        '''
        if string[start]!='(':
            return (None,start)
        branchset,rest=parse_branchset(start+1,end)
        if branchset==None:
            return (None,start)
        if string[rest+1]!=')':
            return (None,start)
        name,rest=parse_name(rest+2,end)
        if name==None:
            return (None,start)
        return ((branchset,name),rest)
    
    def parse_branchset(start,end):
        '''
        BranchSet-> Branch | Branch "," BranchSet
        '''
        branchset=[]
        s = start
        while s<end:
            branch,rest=parse_branch(s,end)
            if branch==None:
                if len(branchset)==0:
                    return (None,Start)
                else:
                    return (branchset,rest)
            else:
                branchset.append(branch)
                if string[rest]!=',':
                    return (branchset,rest)
                else:
                    s=rest+1
        return (None,start)
    
    def parse_branch(start,end):
        '''
        Branch -> Subtree Length
        '''
        subtree,rest=parse_leaf(start,end)
        if subtree==None:        
            return (None,start)
        length,rest = parse_length(rest,end)
        if subtree==None:        
            return (None,start)
        return (subtree,rest)    #FIXME - use length
    
    def parse_name(start,end):
        '''
        Name -> empty | string
        '''
        pos=start
        name=''
        while pos<end:
            if string[pos] in 'dogcat': #FIXME
                name=name+string[pos]
                pos+=1
            else:
                break
        return (name,pos)

    def parse_length(start,end):
        '''
         Length -> empty | ":" number
        '''
        return (1,start)   # FIXME
    
    tree,rest=parse_subtree(0,len(string)-1)
    if tree==None:
        tree,rest=parse_branch(0,len(string)-1)
        if tree==None:
            return (None,0)
    if rest==len(string)-1:
        if string[len(string)-1]==';':
            return (tree,len(string))        
  
    return (None,0)



if __name__=='__main__':
    print (parse('(cat)dog;'))
    print (parse('(dog,cat);'))
    print (parse('dog;'))