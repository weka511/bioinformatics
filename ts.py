#!/usr/bin/env python
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
#    TS Topological  sort

from align import topological_order
from helpers import parse_graph, create_adjacency, format_list

if __name__=='__main__':
     with open(r'C:\Users\Simon\Downloads\rosalind_ts.txt') as f:
          g = parse_graph(f)

          _,_,adj = create_adjacency(g,back=False)
          for k,v in adj.items():
               if k in v:
                    v.remove(k)

          t = topological_order(adj)
          print (format_list(t))
