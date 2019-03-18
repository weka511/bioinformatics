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
#    gs 	General Sink

from graphs import scc,dfs

def gs(edges):
       n,_            = edges[0]
       ccnum,adj,adjr = scc(edges)
       for component in ccnum:
              visited = dfs(adj,n,sequence=(i for i in [component]))
              if len(visited)==n:
                     return component
       return -1

if __name__=='__main__':
       graphs = [
           [(3, 2),
            (3 ,2),
            (2, 1)],
          
           [(3, 2),
            (3, 2),
            (1, 2)]   
       ]
       print (' '.join([str(gs(edges)) for edges in graphs]))