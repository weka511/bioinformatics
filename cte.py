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
    return -1

if __name__=='__main__':
    k = 2
    gs = [[[4, 5],
           [2, 4, 2],
           [3, 2, 1],
           [1, 4, 3],
           [2, 1, 10],
           [1, 3, 4]],

          [[4, 5],
           [3, 2, 1],
           [2, 4, 2],
           [4, 1, 3],
           [2, 1, 10],
           [1, 3, 4]]]
    
    print (' '.join([str(cte(g)) for g in gs]))