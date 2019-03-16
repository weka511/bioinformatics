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

#    NWC Negative Weight Cycle
#

from graphs import bf
        

if __name__=='__main__':
  from helpers import create_list    

  #print(' '.join([str(i) for i in bf(create_list(ext=1))]))
  
  #edges = [(9, 13),
           #(1, 2, 10),
           #(3, 2, 1),
           #(3, 4, 1),
           #(4, 5, 3),
           #(5, 6, -1),
           #(7, 6, -1),
           #(8, 7, 1),
           #(1, 8, 8),
           #(7, 2, -4),
           #(2, 6, 2),
           #(6, 3, -2),
           #(9, 5 ,-10),
           #(9, 4, 7)]
  
  #print (' '.join(str(x) for x in bf(edges) ))
