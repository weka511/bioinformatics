# -*- coding: utf-8 -*-
# snarfed from https://github.com/mutux/Ukkonen-s-Suffix-Tree-Algorithm/blob/master/suffixtree.py
class Node:

    __num__ = -1

    def __init__(self, parentkey, outedges, suffixlink=None):
        self.parentkey = parentkey
        self.outedges = outedges
        self.suffixlink = suffixlink
        Node.__num__ += 1
        self.id = Node.__num__

    def getoutedges(self):
        return self.outedges

    def setoutedge(self, key, param):
        (anode, label_start_index, label_end_index, bnode)=param
        if self.outedges is None:
            self.outedges = {}
        self.outedges[key] = (anode, label_start_index, label_end_index, bnode)

    def getoutedge(self, key):
        if key in self.outedges:
            return self.getoutedges()[key]
        else:
            return None

    def getparenkey(self):
        return self.parentkey

    def setparentkey(self, parentkey):
        self.parentkey = parentkey

    def getsuffixlink(self):
        return self.suffixlink

    def setsuffixlink(self, node):
        self.suffixlink = node

    def getid(self):
        return self.id

    @staticmethod
    def __count__(rnode, chars, v, ed='#'):
        total = 0
        l = len(chars)
        if rnode.getoutedges()==None: return 0
        edges = rnode.getoutedges().items()
        #print (len(edges))
        for key,value in edges:
            _,start,end,linked = value
            fin = end+1 if type(end)==int else -1
            suffix = chars[start:fin ]
            #print (start,end, suffix)
            total+= len(suffix)
            total += Node.__count__(linked,chars,v,ed=ed)
        return total
        
    
    @staticmethod
    def __draw__(rnode, chars, v, ed='#'):
        l = len(chars)
        edges = rnode.getoutedges().items()
        nogc = []
        hasgc = []
        gc = []
        maxlen = len(chars) + 6
        for edg in edges:
            if v == 0:
                if edg[1][3].getoutedges() is None:
                    nogc.append(edg)
                else:
                    hasgc.append(edg)
            else:
                if edg[1][3].getoutedges() is None:
                    hasgc.append(edg)
                else:
                    nogc.append(edg)
        gc.extend(hasgc)
        gc.extend(nogc)
        for k, (parent, s, t, node) in gc:
            if ed == '#':
                if t == '#':
                    t = l
            else:
                if t == '#':
                    t = ed
            linkid = ''
            if node.getsuffixlink() is not None:
                linkid = '->' + str(node.getsuffixlink().getid())

            if v == 0:
                print (" " * maxlen * v + '|')
                print (" " * maxlen * v + '|' + ' ' * 3 + chars[s:t + 1])
                print ('+' + " " * maxlen * v + '-' + '-' * (maxlen - 1) + '● ' + '(' + str(node.getid()) + linkid + ')')
            else:
                print ('|' + " " * maxlen * v + '|')
                print ('|' + " " * maxlen * v + '|' + ' ' * 3 + chars[s:t + 1])
                print ('|' + " " * maxlen * v + '+' + '-' * (maxlen - 1) + '● ' + '(' + str(node.getid()) + linkid + ')')
            if node.getoutedges() is not None:
                Node.__draw__(node, chars, v + 1, ed)

    @staticmethod
    def draw(root, chars, ed='#'):
        #print ('\n', chars, '\n● (0)')
        v = 0
        Node.__draw__(root, chars, v, ed)

    @staticmethod
    def count(root, chars, ed='#'):
        #print ('\n', chars, '\n● (0)')
        v = 0
        return Node.__count__(root, chars, v, ed)
        
def build(chars, regularize=False):
    root = Node(None, None, None)
    actnode = root
    actkey = ''
    actlen = 0
    remainder = 0  # used for splitting
    ind = 0
    while ind < len(chars):
        ch = chars[ind]
        if remainder == 0:
            if actnode.getoutedges() is not None and ch in actnode.getoutedges():
                actkey = ch
                actlen = 1
                remainder = 1
                anode, start, end, bnode = actnode.getoutedge(actkey)
                if end == '#':
                    end = ind
                if end - start + 1 == actlen:
                    actnode = actnode.getoutedge(actkey)[3]
                    actkey = ''
                    actlen = 0
            else:
                aleaf = Node(None, None, None)
                aedge = (actnode, ind, '#', aleaf)
                aleaf.setparentkey((actnode, chars[ind]))
                actnode.setoutedge(chars[ind], aedge)
        else:
            if actkey == '' and actlen == 0:  # compare on node
                if ch in actnode.getoutedges():
                    actkey = ch
                    actlen = 1
                    remainder += 1
                else:
                    remainder += 1
                    remainder, actnode, actkey, actlen = unfold(root, chars, ind, remainder, actnode, actkey, actlen)
            else:  # compare on edge
                anode, start, end, bnode = actnode.getoutedge(actkey)
                if end == '#':
                    end = ind
                compareposition = start + actlen
                if chars[compareposition] != ch:
                    remainder += 1
                    remainder, actnode, actkey, actlen = unfold(root, chars, ind, remainder, actnode, actkey, actlen)
                else:
                    if compareposition < end:  # on edge
                        actlen += 1
                        remainder += 1
                    else:  # on node
                        remainder += 1
                        actnode = actnode.getoutedge(actkey)[3]
                        if compareposition == end:
                            actlen = 0
                            actkey = ''
                        else:
                            actlen = 1
                            actkey = ch
        ind += 1
        if ind == len(chars) and remainder > 0:
            if regularize:
                chars = chars + '$'
    return root, chars


def unfold(root, chars, ind, remainder, actnode, actkey, actlen):
    prenode = None
    while remainder > 0:
        remains = chars[ind - remainder + 1:ind + 1]
        actlen_re = len(remains) - 1 - actlen
        actnode, actkey, actlen, actlen_re = hop(ind, actnode, actkey, actlen, remains, actlen_re)
        lost, actnode, actkey, actlen, actlen_re = step(chars, ind, actnode, actkey, actlen, remains, actlen_re)
        if lost:
            if actlen == 1 and prenode is not None and actnode is not root:
                prenode.setsuffixlink(actnode)
            return remainder, actnode, actkey, actlen
        if actlen == 0:
            if remains[actlen_re] not in actnode.getoutedges():
                aleaf = Node(None, None, None)
                aedge = (actnode, ind, '#', aleaf)
                aleaf.setparentkey((actnode, chars[ind]))
                actnode.setoutedge(chars[ind], aedge)
        else:  # on edge
            anode, start, end, bnode = actnode.getoutedge(actkey)
            if remains[actlen_re + actlen] != chars[start + actlen]:
                # split
                anode, start, end, bnode = actnode.getoutedge(actkey)
                newnode = Node(None, None, None)
                halfedge1 = (actnode, start, start + actlen - 1, newnode)
                halfedge2 = (newnode, start + actlen, end, bnode)
                actnode.setoutedge(actkey, halfedge1)
                newnode.setparentkey((actnode, actkey))
                newnode.setoutedge(chars[start + actlen], halfedge2)
                aleaf = Node(None, None, None)
                aedge = (newnode, ind, '#', aleaf)
                aleaf.setparentkey((newnode, chars[ind]))
                newnode.setoutedge(chars[ind], aedge)
            else:
                return remainder, actnode, actkey, actlen
        if prenode is not None and 'aleaf' in locals() and aleaf.getparenkey()[0] is not root:
            prenode.setsuffixlink(aleaf.getparenkey()[0])
        if 'aleaf' in locals() and aleaf.getparenkey()[0] is not root:
            prenode = aleaf.getparenkey()[0]
        if actnode == root and remainder > 1:
            actkey = remains[1]
            actlen -= 1
        if actnode.getsuffixlink() is not None:
            actnode = actnode.getsuffixlink()
        else:
            actnode = root
        remainder -= 1
    return remainder, actnode, actkey, actlen


def step(chars, ind, actnode, actkey, actlen, remains, ind_remainder):
    rem_label = remains[ind_remainder:]
    if actlen > 0:
        anode, start, end, bnode = actnode.getoutedge(actkey)
        if end == '#':
            end = ind
        edgelabel = chars[start:end + 1]
        if edgelabel.startswith(rem_label):
            actlen = len(rem_label)
            actkey = rem_label[0]
            return True, actnode, actkey, actlen, ind_remainder
    else:
        # on node
        if ind_remainder < len(remains) and remains[ind_remainder] in actnode.getoutedges():
            anode, start, end, bnode = actnode.getoutedge(remains[ind_remainder])
            if end == '#':
                end = ind
            edgelabel = chars[start:end + 1]
            if edgelabel.startswith(rem_label):
                actlen = len(rem_label)
                actkey = rem_label[0]
                return True, actnode, actkey, actlen, ind_remainder
    return False, actnode, actkey, actlen, ind_remainder


def hop(ind, actnode, actkey, actlen, remains, ind_remainder):
    if actlen == 0 or actkey == '':
        return actnode, actkey, actlen, ind_remainder
    anode, start, end, bnode = actnode.getoutedge(actkey)
    if end == '#':
        end = ind
    edgelength = end - start + 1
    while actlen > edgelength:
        actnode = actnode.getoutedge(actkey)[3]
        ind_remainder += edgelength
        actkey = remains[ind_remainder]
        actlen -= edgelength
        anode, start, end, bnode = actnode.getoutedge(actkey)
        if end == '#':
            end = ind
        edgelength = end - start + 1
    if actlen == edgelength:
        actnode = actnode.getoutedge(actkey)[3]
        actkey = ''
        actlen = 0
        ind_remainder += edgelength
    return actnode, actkey, actlen, ind_remainder


if __name__ == "__main__":
    docs = ['abcabxabcd', 'dedododeeodoeodooedeeododooodoede$', 'ooooooooo', 'mississippi']
    for text in docs:
        tree, pst = build(text, regularize=True)
        Node.draw(tree, pst, ed='#')