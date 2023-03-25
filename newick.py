#!/usr/bin/env python
#  Copyright (C) 2017-2023 Greenweaves Software Limited

#  This is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this software.  If not, see <http://www.gnu.org/licenses/>

#   Grammar from https://en.wikipedia.org/wiki/Newick_format
#      Tree -> Subtree ";" | Branch ";"
#      Subtree -> Leaf | Internal
#      Leaf -> Name
#      Internal -> "(" BranchSet ")" Name
#      BranchSet-> Branch | Branch "," BranchSet
#      Branch -> Subtree Length
#      Name -> empty | string
#      Length -> empty | ":" number

from string import ascii_letters, digits, whitespace
from re     import findall


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
        self.parent = None

        Node.next_id +=1

    def add(self,node):
        #print ('Added {0}({1}) to {2}({3})'.format(node.id,node.level,self.id,self.level))
        self.nodes.append(node)
        node.parent = self

    def iterate(self,path=[]):
        yield self,path
        new_path=path[:]
        new_path.append(self.id)
        for node in self.nodes:
            for rn in node.iterate(new_path):
                yield rn

class Tokenizer:
    UNDEFINED         = -2
    WHITESPACE        = -1
    SEMICOLON         = 0
    OPEN_PARENTHESIS  = 1
    CLOSE_PARENTHESIS = 2
    COMMA             = 3
    COLON             = 4
    NAME              = 5
    NUMBER            = 6
    def __init__(self):
        self.symbols={
            ';' : Tokenizer.SEMICOLON,
            '(' : Tokenizer.OPEN_PARENTHESIS,
            ')' : Tokenizer.CLOSE_PARENTHESIS,
            ',' : Tokenizer.COMMA,
            ':' : Tokenizer.COLON,
            '_' : Tokenizer.NAME,
            '.' : Tokenizer.NUMBER,
            }
        for i in ascii_letters:
            self.symbols[i] = Tokenizer.NAME
        for i in digits:
            self.symbols[i] = Tokenizer.NUMBER
        for i in whitespace:
            self.symbols[i] = Tokenizer.WHITESPACE

    def get_token(self,ch):
        try:
            return self.symbols[ch]
        except KeyError as e:
            print('Invalid character: {0} {1}'.format(ch,e))
            return Tokenizer.UNDEFINED

    def tokenize(self,text):
        def long_token(target,pos):
            start = pos-1
            token = target
            while token==target:
                ch=text[pos]
                token = self.get_token(ch)
                pos+=1
            return start,pos
        pos = 0
        end = len(text)
        while pos<end:
            ch=text[pos]
            token = self.get_token(ch)
            pos+=1
            if token==Tokenizer.UNDEFINED:
                continue
            if token==Tokenizer.WHITESPACE:
                continue
            elif token<Tokenizer.NAME:
                yield token,pos-1,pos,None
            elif token==Tokenizer.NAME:
                start,pos=long_token(token,pos)
                pos-=1
                yield token,start,pos,text[start:pos]
            elif token==Tokenizer.NUMBER:
                start,pos=long_token(token,pos)
                pos-=1
                yield token,start,pos,self.value(text[start:pos])
    def value(self,text):
        try:
            return int(text)
        except ValueError:
            try:
                return float(text)
            except ValueError as e:
                print ('Could not parse length: ',e)

class Parser():

    def __init__(self,tokenizer):
        self.tokenizer=tokenizer

    def parse(self,text):
        ended   = False
        stack   = []
        current = None
        tree    = None
        lookup  = {}
        Node.next_id = 0

        for token,start,pos,value in self.tokenizer.tokenize(text):
            if ended:
                print ('Characters beyond end',token,start,pos,text[staself.symbolsrt:pos])
            #else:
                #print (token,start,pos,text[start:pos])
            if token==Tokenizer.SEMICOLON:
                ended = True
            elif token==Tokenizer.OPEN_PARENTHESIS:
                branchset = Node(len(stack),start)
                if len(stack)>0:
                    stack[len(stack)-1].add(branchset)
                stack.append(branchset)
                if tree==None:
                    tree = stack[0]
                current=Node(len(stack),start)       # first branch
                stack[len(stack)-1].add(current)
            elif token==Tokenizer.CLOSE_PARENTHESIS:
                current=stack.pop()
                current.end = pos
            elif token==Tokenizer.COMMA:
                current=Node(len(stack),start)
                stack[len(stack)-1].add(current)
            elif token==Tokenizer.COLON:
                pass
            elif token==Tokenizer.NUMBER:
                current.length=value
            elif token==Tokenizer.NAME:
                current.name=value
                lookup[value]=current
        return (tree,lookup)

# newick_to_adjacency_list

def newick_to_adjacency_list(T,return_root=False):
    Adj       = {}
    Stack     = [[]]
    top       = []
    tokenizer = Tokenizer()
    root      = None
    for token,start,pos,value in tokenizer.tokenize(T):
        if token==Tokenizer.OPEN_PARENTHESIS:
            Stack.append([])
        elif token==Tokenizer.NAME:
            Stack[-1].append(value)
            Adj[value] = top
            top        = []
            root       = value
        elif token==Tokenizer.CLOSE_PARENTHESIS:
            top = Stack.pop()
        elif token==Tokenizer.SEMICOLON:
            pass

    if return_root:
        return Adj,root
    else:
        return Adj

class Hierarchy:

    '''This class parses a Newick string into a tree'''

    def __init__(self,newick,start=0):
        '''
        parse

        Parse a new string into a hierarchy

        # snarfed from https://stackoverflow.com/questions/51373300/how-to-convert-newick-tree-format-to-a-tree-like-hierarchical-object
        '''
        tokens = findall(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

        def recurse(nextid = start, parentid = -1): # one node
            thisid = nextid;
            children = []

            name, length, delim, ch = tokens.pop(0)
            if ch == "(":
                while ch in "(,":
                    node, ch, nextid = recurse(nextid+1, thisid)
                    children.append(node)
                name, length, delim, ch = tokens.pop(0)
            return {"id": thisid, "name": name, "length": float(length) if length else None,
                    "parentid": parentid, "children": children}, delim, nextid

        self.tree = recurse()[0]

    def create_adj(self):
        adj = {}
        def dfs(tree):
            id       = tree['id']
            name     = tree['name']
            children = tree['children']
            parentid = tree['parentid']
            if len(name)==0:
                adj[id]=[]
            if parentid>-1:
                adj[parentid].append(id if len(name)==0 else name)
            for child in children:
                dfs(child)
        dfs(self.tree)
        return adj


if __name__=='__main__':
    def display(p):
        for (node,path) in p.iterate():
            spaces = ''.join(['-' for i in range(node.level)])
            print ('{5} {0}ID={1}, Name={2}, level={3}, length={4}'.format(spaces,node.id, node.name,node.level,node.length,path))

    parser = Parser(Tokenizer())

    s = '''
        ((raccoon:19.1.9959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,
        ((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);
        '''
    p,lookup=parser.parse(s)
    print ('----------')
    display(p)
    print (lookup)
    #display (parser.parse('(cat)dog;'))
    #print (parse('dog;'))
    #display (parser.parse('(dog,cat);'))
    #pp=parser.parse(
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
    #display(pp)
    #for node in pp.iterate():
        #spaces = ''.join(['-' for i in range(node.level)])
        #print ('{0}{1} {2} {3}'.format(spaces,node.id, node.name,node.level))
