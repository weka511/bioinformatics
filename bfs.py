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
    for (a,b) in links:
        if a in result:
            result[a].append(b)
        else:
            result[a]=[b]
        
    return result


def BestFirstSearch(tree,start,path=[],sub_total=0):
    result=[(start,sub_total)]
    if start in tree:
        for node in tree[start]:
            if not node in path:
                result=result+BestFirstSearch(tree,node,sub_total=sub_total+1,path=path+[node])
    if len(path)==0:
        n=max(a for (a,b) in result)
        result.sort()
        final=[]
        i=0
        b,c = result[i]
        print (result)
        for a in range(1,n):          
            if a<b:
                final.append((a,-1))
            else:
                i+=1
                b0,c0 = result[i]
                while b0==b:
                    if c0<c:
                        c=c0
                    i+=1
                    b0,c0 = result[i]
                final.append((b,c))
                b=b0
                c=c0
                if i>len(result): break
        return [b for (a,b) in final+[(a+1,c)]]
    else:
        return result

if __name__=='__main__':
    print(
        BestFirstSearch(tree([
            (6, 6),
            (4, 6),
            (6, 5),
            (4, 3),
            (3, 5),
            (2, 1),
            (1, 4)]),1))

    