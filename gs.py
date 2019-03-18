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

from graphs import gs



if __name__=='__main__':
       
       from helpers import create_strings
       
       #graphs = [
           #[(3, 2),
            #(3 ,2),
            #(2, 1)],
          
           #[(3, 2),
            #(3, 2),
            #(1, 2)]   
       #]
       
       graphs = []
       edges = []
       for s in create_strings(ext=2):
              if len(s)==0:
                     if len(edges)>0:
                            graphs.append(edges)
                            edges=[]
              else:
                     values = [int(x) for x in s.split(" ")]
                     if len(values)>1:
                            edges.append(values)
       if len(edges)>0:
              graphs.append(edges)     
              
       print (' '.join([str(gs(edges)) for edges in graphs]))