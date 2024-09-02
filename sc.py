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
import os
from time import time
from graphs import sc
from helpers import create_list

def generate_graphs():
   state = 0
   for entry in create_list(path='./data'):
      match(state):
         case 0:
            print (f'There are {entry[0]} graphs')
            state = 1
         case 1:
            edges = [entry]
            state = 2
            _,n = entry
         case 2:
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
      for graph in generate_graphs():
         print (sc(graph))


      # print (' '.join([str(sc(g)) for g in graphs]))

   elapsed = time() - start
   minutes = int(elapsed/60)
   seconds = elapsed - 60*minutes
   print (f'Elapsed Time {minutes} m {seconds:.2f} s')
