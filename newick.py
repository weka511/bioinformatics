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
    def iterate(self):
        yield self
        for node in self.nodes:
            for rn in node.iterate():
                yield rn

    
def parse(newick_text):   
    def find_parens(start,end,tree,stack):
        for i in range(start,end):
            if string[i]=='(':
                top = stack[len(stack)-1]               
                node=Node(top.level+1)
                node.start=i
                top.add(node)
                stack.append(node)
            elif string[i]==')':
                top = stack.pop()
                top.end=i
            
        return tree
    
    string= ''.join(newick_text.split()) 
    print (string)
    tree = Node(0)
    
    return find_parens(0,len(string),tree,[tree])
    

if __name__=='__main__':
    s = '''
        ((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,
        ((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);
        '''
    t=string= ''.join(s.split())
    for node in parse(s).iterate():
            print (node.level, node.start, node.end,t[node.start:node.end+1])
    #print (parse('(cat)dog;'))
    #print (parse('dog;'))
    #print (parse('(dog,cat);'))
    #pp=parse(
        #'''
        #(((Abantias_ibis,Phormictopus_fissipes),(((Alauda_guineti,Certhia_geniculata),(Balaena_paradoxus,Mustela_riparia)),
        #((Circus_guentheri,Oligodon_cyanus),(Gazella_ruthveni,Terpsihone_dentatus)))),
        #((((Aegialifes_duplex,(Athene_sinensis,Underwoodisaurus_subglobosa)),(Damon_cinclus,Nemachilus_sepsoides)),
        #((Pelusios_trigonopodus,Porzana_minutus),Uroplatus_gallicus)),Tetrao_cristatus),
        #((((((((((((((((((((((Ahaetulla_nasuta,Cypselus_drapiezii),Rosalia_hypomelus),
        #(Gallinago_scutulata,Sturnus_conicus)),Synthliboramphus_peregusna),(((((Alectoris_tridactylum,
        #(((Alloporus_torquatus,((((Argynnis_flavolineata,Phrynops_flavigularis),
        #(((Chelydra_blythi,(Lepidobatrachus_pardus,Pelecanus_lutris)),(Ortigometra_colchicus,Rhacodactylus_dolosus)),
        #Pseudorca_maculata)),Gazella_colombianus),Physignathus_noctua)),(Eucratoscelus_taxus,Nemachilus_chuatsi)),Burhinus_lutra)),
        #Saga_chinensis),Hadogenes_verrucosus),Eryx_carinata),(Anodonta_zagrosensis,((Leiolepis_perdix,Lepus_zagrosensis),
        #Nyroca_trigonopodus)))),(Gekko_aeruginosus,Mochlus_Bernicla)),Castor_leucomelas),Rhesus_caninus),
        #(Gyps_ovata,((Hirundo_rapax,Xenophrys_fissipes),Physignathus_collaris))),(Bronchocela_chrysaetus,Ortigometra_jubata)),
        #(Cyriopagopus_elaphus,Middendorffinaia_plumipes)),Chelus_lepturus),Kassina_mexicanum),((Chalcides_dendrophila,Leiocephalus_unicolor),
        #Ovis_musculus)),Saxicola_nasuta),Phasianus_modestus),Chettussia_pholeter),Heterodon_physalus),((Mustela_teguixin,Trachemys_marmoratus),
        #Scaphiopus_nelsonii)),Certhia_flavirufa),(Aphonopelma_gratiosa,Rhodostethia_macrops)),(((Circaetus_olivacea,
        #((((Coturnix_pyromelana,Porzana_corticale),Hyla_ussuriensis),(Falcipennis_taxispilota,(Notophthalmus_macularius,Rissa_unicus))),
        #Upupa_himantopus)),Circus_homeana),Riparia_ferruginea)));
        #''')
    #for node in pp.iterate():
        #print (node.level, node.start, node.end)