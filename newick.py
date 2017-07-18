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

import string


class Node:
    next_id = 0
    def __init__(self,level,start):
        self.start  = start
        self.end    = 0
        self.nodes  = []
        self.level  = level
        self.length = 1
        self.name   = '-'
        self.id     = Node.next_id
        Node.next_id +=1
        print ("Node: ",self.id, self.level)
        
    def add(self,node):
        print ('Added {0}({1}) to {2}({3})'.format(node.id,node.level,self.id,self.level))
        self.nodes.append(node)

    def iterate(self):
        yield self
        for node in self.nodes:
            #yield node
            for rn in node.iterate():
                yield rn

WHITESPACE        = -1
SEMICOLON         = 0
OPEN_PARENTHESIS  = 1
CLOSE_PARENTHESIS = 2
COMMA             = 3
COLON             = 4
NAME              = 5
NUMBER            = 6  

class Parser():
    def __init__(self):
        self.symbols={
            ';' : SEMICOLON,
            '(' : OPEN_PARENTHESIS,
            ')' : CLOSE_PARENTHESIS,
            ',' : COMMA,
            ':' : COLON,
            '_' : NAME,
            '.' : NUMBER,
            }
        for i in string.ascii_letters:
            self.symbols[i] = NAME
        for i in string.digits:
            self.symbols[i] = NUMBER
        for i in string.whitespace:
            self.symbols[i] = WHITESPACE

    def tokenize(self,text):        
        def long_token(target,pos):
            start = pos-1
            token = target
            while token==target:
                ch=text[pos]
                token = self.symbols[ch]
                pos+=1
            return start,pos
        pos = 0
        end = len(text)
        while pos<end:
            ch=text[pos]
            token = self.symbols[ch]
            pos+=1
            if token==WHITESPACE:
                continue
            elif token<NAME:
                yield token,pos-1,pos
            elif token==NAME or token==NUMBER:
                start,pos=long_token(token,pos)
                pos-=1
                yield token,start,pos
                        
    def parse(self,text):        
        ended = False
        stack = []
        current = None
        tree = None
        for token,start,pos in self.tokenize(text):
            if ended:
                print ('Characters beyond end',token,start,pos,text[staself.symbolsrt:pos])
            else:
                print (token,start,pos,text[start:pos])
            if token==SEMICOLON:
                print (len(stack))
                ended = True
            elif token==OPEN_PARENTHESIS:
                new_paren = Node(len(stack),start)
                if len(stack)>0:
                    stack[len(stack)-1].add(new_paren)
                stack.append(new_paren) # the branchset
                if tree==None:
                    tree = stack[0]
                current=Node(len(stack),start)       # first branch
                stack[len(stack)-1].add(current)
            elif token==CLOSE_PARENTHESIS:
                current=stack.pop()
                current.end = pos
            elif token==COMMA:
                current=Node(len(stack),start)   
                stack[len(stack)-1].add(current)
            elif token==COLON:
                pass
            elif token==NUMBER:
                pass
            elif token==NAME:
                current.name=text[start:pos]
        return tree


if __name__=='__main__':
    parser = Parser()
    
    s = '''
        ((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,
        ((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);
        '''
    p=parser.parse(s)
    print ('----------')
    for node in p.iterate():
        spaces = ''.join(['-' for i in range(node.level)])
        print ('{0}{1} {2} {3}'.format(spaces,node.id, node.name,node.level))
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