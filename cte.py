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

# CTE Shortest Cycle Through a Given Edge

def cte(g):
    open_edges     = g[1:]
    a,b,w          = open_edges[0]
    open_edges     = open_edges[1:]
    partial_cycles = [(b,w)]
    cycles         = []
    start          = a
    for i in range(len(open_edges)):
        new_cycles   = []
        closed_edges = {}
        for a,b,w in open_edges:
            for end,n in partial_cycles:
                if end==a:
                    new_cycles.append((b,w+n))
        partial_cycles=[]
        for b,w in new_cycles:
            if b==start:
                cycles.append((b,w))
            else:
                partial_cycles.append((b,w))
        if len(partial_cycles)==0:
            return min([w for (_,w) in cycles]) if len(cycles)>0 else -1
    return -1

if __name__=='__main__':
    k = 2
    gs = [[(4, 5),
           (2, 4, 2),
           (3, 2, 1),
           (1, 4, 3),
           (2, 1, 10),
           (1, 3, 4)],

          [(4, 5),
           (3, 2, 1),
           (2, 4, 2),
           (4, 1, 3),
           (2, 1, 10),
           (1, 3, 4)]]
    
    print (' '.join([str(cte(g)) for g in gs]))