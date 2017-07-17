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

#def parse(string):
    #'''
    #Tree -> Subtree ";" | Branch ";
    #'''
    #def parse_subtree(start,end):
        #'''
         #Subtree -> Leaf | Internal
        #'''
        #print ('Parse subtree {0}[{1}:{2}]'.format(string,start,end))
        #subtree,rest=parse_leaf(start,end)
        #if subtree==None:
            #subtree,rest=parse_internal(start,end)
        #if subtree==None:
            #return (None,start)
        #print ('Subtree: {0},pos={1}'.format(subtree,rest))
        #return (subtree,rest)
    
    #def parse_leaf(start,end):
        #'''
         #Leaf -> Name
         #'''
        #print ('Parse leaf {0}[{1}:{2}]'.format(string,start,end))
        #leaf,rest=parse_name(start,end)
        #if leaf==None:
            #return (None,start)
        #print ('Leaf: {0},pos={1}'.format(leaf,rest))
        #return (leaf,rest)
    
    #def parse_internal(start,end):
        #'''
        #Internal -> "(" BranchSet ")" Name
        #'''
        #print ('Parse internal {0}[{1}:{2}]'.format(string,start,end))
        #if string[start]!='(':
            #return (None,start)
        #branchset,rest=parse_branchset(start+1,end)
        #if branchset==None:
            #return (None,start)
        #if string[rest+1]!=')':
            #return (None,start)
        #name,rest=parse_name(rest+2,end)
        #if name==None:
            #return (None,start)
        #print ('Internal: {0},pos={1}'.format(branchset,rest))
        #return ((branchset,name),rest)
    
    #def parse_branchset(start,end):
        #'''
        #BranchSet-> Branch | Branch "," BranchSet
        #'''
        #print ('Parse branchset {0}[{1}:{2}]'.format(string,start,end))
        #branchset=[]
        #s = start
        #while s<end:
            #branch,rest=parse_branch(s,end)
            #if branch==None:
                #if len(branchset)==0:
                    #return (None,Start)
                #else:
                    #print ('Branchset: {0},pos={1}'.format(branchset,rest))
                    #return (branchset,rest)
            #else:
                #branchset.append(branch)
                #if string[rest]!=',':
                    #print ('Branchset: {0},pos={1}'.format(branchset,rest))
                    #return (branchset,rest)
                #else:
                    #s=rest+1
        #return (None,start)
    
    #def parse_branch(start,end):
        #'''
        #Branch -> Subtree Length
        #'''
        #print ('Parse branch {0}[{1}:{2}]'.format(string,start,end))
        #subtree,rest=parse_leaf(start,end)
        #if subtree==None:        
            #return (None,start)
        #length,rest = parse_length(rest,end)
        #if subtree==None:        
            #return (None,start)
        #print ('Branch: {0},pos={1}'.format(branch,rest))
        #return (subtree,rest)    #FIXME - use length
    
    #def parse_name(start,end):
        #'''
        #Name -> empty | string
        #'''
        #print ('Parse name {0}[{1}:{2}]'.format(string,start,end))
        #pos=start
        #name=''
        #while pos<end:
            #if string[pos] in 'dogcat': #FIXME
                #name=name+string[pos]
                #pos+=1
            #else:
                #break
            
        #print ('Name: {0},pos={1}'.format(name,pos))
        #return (name,pos)

    #def parse_length(start,end):
        #'''
         #Length -> empty | ":" number
        #'''
        #print ('Parse length {0}[{1}:{2}]'.format(string,start,end))
        #return (1,start)   # FIXME
    
    #print ('Parse {0}'.format(string))
    #tree,rest=parse_subtree(0,len(string)-1)
    #if tree==None:
        #tree,rest=parse_branch(0,len(string)-1)
        #if tree==None:
            #return (None,0)
    #if rest==len(string)-1:
        #if string[len(string)-1]==';':
            #return (tree,len(string))        
  
    #return (None,0)

class Node:
    def __init__(self,level):
        self.start = 0
        self.end   = 0
        self.nodes = []
        self.level = level
        
    def add(self,node):
        self.nodes.append(node)
    def __str__(self):
        children=''
        for node in self.nodes:
            children = children + str(node)
        header=''
        for i in range(self.level):
            header=header+'-'
        return '{0}{1},{2}\n{3}'.format(header,self.start,self.end,children) 
    
def parse(string):
    
    def find_parens(start,end,tree,stack):
        for i in range(start,end):
            if string[i]=='(':
                top = stack[len(stack)-1]               
                node=Node(top.level+1)
                node.start=i
                top.add(node)
                stack.append(node)
            elif string[i]==')':
                stack.pop().end=i
            
        return tree
    
    print (string)
    tree = Node(0)
    
    return find_parens(0,len(string),tree,[tree])
    

if __name__=='__main__':
    print (parse('((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);'))
    #print (parse('(cat)dog;'))
    #print (parse('dog;'))
    #print (parse('(dog,cat);'))
    print (parse(
        '''
        (((Abantias_ibis,Phormictopus_fissipes),(((Alauda_guineti,Certhia_geniculata),(Balaena_paradoxus,Mustela_riparia)),
        ((Circus_guentheri,Oligodon_cyanus),(Gazella_ruthveni,Terpsihone_dentatus)))),
        ((((Aegialifes_duplex,(Athene_sinensis,Underwoodisaurus_subglobosa)),(Damon_cinclus,Nemachilus_sepsoides)),
        ((Pelusios_trigonopodus,Porzana_minutus),Uroplatus_gallicus)),Tetrao_cristatus),
        ((((((((((((((((((((((Ahaetulla_nasuta,Cypselus_drapiezii),Rosalia_hypomelus),
        (Gallinago_scutulata,Sturnus_conicus)),Synthliboramphus_peregusna),(((((Alectoris_tridactylum,
        (((Alloporus_torquatus,((((Argynnis_flavolineata,Phrynops_flavigularis),(((Chelydra_blythi,(Lepidobatrachus_pardus,Pelecanus_lutris)),(Ortigometra_colchicus,Rhacodactylus_dolosus)),Pseudorca_maculata)),Gazella_colombianus),Physignathus_noctua)),(Eucratoscelus_taxus,Nemachilus_chuatsi)),Burhinus_lutra)),Saga_chinensis),Hadogenes_verrucosus),Eryx_carinata),(Anodonta_zagrosensis,((Leiolepis_perdix,Lepus_zagrosensis),Nyroca_trigonopodus)))),(Gekko_aeruginosus,Mochlus_Bernicla)),Castor_leucomelas),Rhesus_caninus),(Gyps_ovata,((Hirundo_rapax,Xenophrys_fissipes),Physignathus_collaris))),(Bronchocela_chrysaetus,Ortigometra_jubata)),(Cyriopagopus_elaphus,Middendorffinaia_plumipes)),Chelus_lepturus),Kassina_mexicanum),((Chalcides_dendrophila,Leiocephalus_unicolor),Ovis_musculus)),Saxicola_nasuta),Phasianus_modestus),Chettussia_pholeter),Heterodon_physalus),((Mustela_teguixin,Trachemys_marmoratus),Scaphiopus_nelsonii)),Certhia_flavirufa),(Aphonopelma_gratiosa,Rhodostethia_macrops)),(((Circaetus_olivacea,((((Coturnix_pyromelana,Porzana_corticale),Hyla_ussuriensis),(Falcipennis_taxispilota,(Notophthalmus_macularius,Rissa_unicus))),Upupa_himantopus)),Circus_homeana),Riparia_ferruginea)));
Certhia_geniculata Eryx_carinata
'''))