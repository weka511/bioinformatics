#    Copyright (C) 2017 Greenweaves Software Pty Ltd
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
#    bfs 	Beadth First search

def tree(links):
    result={}
    max_node=-1
    for (a,b) in links:
        if a>max_node: max_node=a
        if b>max_node: max_node=b
        
        if a in result:
            result[a].append(b)
        else:
            result[a]=[b]
        
    return (max_node,result)

def ShortestDistances(tree):
    result=[]
    def BreadthFirstSearch(start,dist=0):
        result[start-1]=dist
        dist+=1
        if start in tree:
            for node in tree[start]:
                if result[node-1]<0 or result[node-1]>dist:
                    BreadthFirstSearch(node,dist)   
                    
    max_node,tree=tree
    result=[-1]*max_node
    BreadthFirstSearch(1)
    return result
    


if __name__=='__main__':
    print(
        ShortestDistances(tree([
            #(6, 6),
            (4, 6),
            (6, 5),
            (4, 3),
            (3, 5),
            (2, 1),
            (1, 4)])))

    with open ('c:/Users/Weka/Downloads/rosalind_bfs(4).txt') as f:
        ll=[]
        i=0
        for line in f:
            if i>0:
                xs=line.strip().split()
                ll.append((int(xs[0]),int(xs[1])))
            i+=1
        print (ShortestDistances(tree(ll)))