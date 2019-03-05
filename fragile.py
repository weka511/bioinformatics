#    Copyright (C) 2019 Greenweaves Software Limited
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>
#
# fragile.py # code for Chapter 6 - Are there fragile regions in the human genome?

# BA6B Compute the Number of Breakpoints in a Permutation 

def getBreakPoints(P):
    return len( [(a,b) for (a,b) in zip([0]+P,P+[len(P)+1]) if b!=a+1] )

# BA6C Compute the 2-Break Distance Between a Pair of Genomes

def get_synteny_blocks(a):
    return [j for i in a for j in i] 

def count_synteny_blocks(a):
    return max(abs(b) for b in get_synteny_blocks(a))

# BA6F Implement Chromosome to Cycle

def ChromosomeToCycle(Chromosome):
    Nodes = []
    for i in Chromosome:
        if i> 0:
            Nodes.append(2*i-1)
            Nodes.append(2*i)
        else:
            Nodes.append(-2*i)
            Nodes.append(-2*i-1)
        
    return Nodes

# BA6H Implement ColoredEdges

def ColouredEdges(P):
    Edges = []
    for  Chromosome in P:
        Nodes = ChromosomeToCycle(Chromosome)
        it = iter(Nodes[1:]+[Nodes[0]])
        for i in it:
            Edges.append((i,next(it)))
    return Edges

# BA6I Implement Graph to Genome

def GraphToGenome(GenomeGraph):
    def diff(a,b):
        _,x=a
        y,_=b
        return abs(x-y)
    def build_cycle(pair,dcycles):
        result = [pair]
        while pair in dcycles:
            pair = dcycles[pair]
            result.append(pair)
        return result
    def black_edges(cycle):
        result=[]
        for i in range(len(cycle)):
            a,_ = cycle[i]
            _,b = cycle[i-1]
            result.append((b,a))
        return result
 
    extract = [(a,b,diff(a,b)) for (a,b) in zip([GenomeGraph[-1]]+GenomeGraph[0:],GenomeGraph)]
    gaps    = [(a,b) for (a,b,diff) in extract if diff>1]
    cycles  = [(a,b) for (a,b,diff) in extract if diff==1]
    dcycles = {}
    for (a,b) in cycles:
        dcycles[a]=b
    P       = [build_cycle(pair,dcycles) for _,pair in gaps]
    Q       = [black_edges(p) for p in P]
    next_node = 1
    R = []
    for q in Q:
        r = []
        for a,b in q:
            r.append(next_node if a<b else -next_node)
            next_node+=1
        R.append(r)
    return R