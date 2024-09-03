#!/usr/bin/env python
#    Copyright (C) 2019-2024 Greenweaves Software Limited
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

'''  sc 	Semi-Connected Graph'''

from argparse import ArgumentParser
from time import time
from graphs import sc
from helpers import create_list

def generate_graphs(path='./data'):
   '''
   Parse input file
   '''
   state = 0
   for entry in create_list(path=path):
      match(state):
         case 0:               # Start of file: number of graphs
            state = 1
         case 1:                # Start of graph: number of nodes and adges
            edges = [entry]
            state = 2
            _,n = entry
         case 2:                 # An edge
            edges.append(entry)
            n -= 1
            if n == 0:
               state = 1
               yield edges

if __name__=='__main__':

   start = time()
   parser = ArgumentParser(__doc__)
   parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
   parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
   args = parser.parse_args()
   if args.sample:
      print (sc( [(3, 2),
                  (3 ,2),
                  (2, 1)]))

   if args.rosalind:
      print (' '.join([str(sc(g)) for g in generate_graphs()]))

   elapsed = time() - start
   minutes = int(elapsed/60)
   seconds = elapsed - 60*minutes
   print (f'Elapsed Time {minutes} m {seconds:.2f} s')
