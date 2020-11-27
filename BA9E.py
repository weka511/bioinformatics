#  Copyright (C) 2020 Greenweaves Software Limited
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# BA9E Find the Longest Substring Shared by Two Strings

import argparse
from   helpers import read_strings
from   snp     import SuffixTree
import sys
import time

# See https://en.wikipedia.org/wiki/Longest_common_substring_problem

def LongestSharedSubstring(s,t):
    tree = SuffixTree()
    Leaves = tree.build(s + '$' + t + '#')
    tree.print()

    
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9E Find the Longest Substring Shared by Two Strings')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (LongestSharedSubstring('TCGGTAGATTGCGCCCACTC',
                                 'AGGGGCTCGCAGTGTAAGAA'))
        
    if args.extra:
        #sys.setrecursionlimit(4000)
        Input,Expected  = read_strings('data/LongestSharedSubstring.txt',init=0)
        #print (Input[0])       
        Actual = LongestSharedSubstring(Input[0],Input[1])
        #print (len(Expected[0]),len(Actual))
        #print (Expected[0])
        #print (Actual)        
 
    if args.rosalind:
        pass
    
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')