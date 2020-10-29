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
from helpers import read_strings
from snp import FindLongestRepeat
import sys
import time

#def FindLongestSubstring(text1,text2):
    #tree1 = SuffixTree()
    #tree1.build(text1)
    #edges1 = sorted(list(set(tree1.collectEdges()))) #sort_by_descending_length(tree1.collectEdges()) 
    #for e in edges1:
        #print (e) 
    #print ('-------------')
    #tree2 = SuffixTree()
    #tree2.build(text2)
    #edges2 = sorted(list(set(tree2.collectEdges())))#sort_by_descending_length(tree2.collectEdges())
    #for e in edges2:
        #print (e)
    #best_match=[]
    #iter1 = iter(edges1)
    #iter2 = iter(edges2)
    #edge1 = next(iter1)
    #edge2 = next(iter2)
    #while True:
        #if edge1<edge2:
            #pass
        #elif edge2<edge1:
            #pass
        #else:
            #pass
    
if __name__=='__main__':
    start = time.time()
    parser = argparse.ArgumentParser('BA9E Find the Longest Substring Shared by Two Strings')
    parser.add_argument('--sample',   default=False, action='store_true', help='process sample dataset')
    parser.add_argument('--extra',    default=False, action='store_true', help='process extra dataset')
    parser.add_argument('--rosalind', default=False, action='store_true', help='process Rosalind dataset')
    args = parser.parse_args()
    if args.sample:
        print (FindLongestRepeat('TCGGTAGATTGCGCCCACTC',
                                 'AGGGGCTCGCAGTGTAAGAA'))
        
    if args.extra:
        sys.setrecursionlimit(2000)
        Input,Expected  = read_strings('data/LongestSharedSubstring.txt',init=0)
        #print (Input[0])       
        Actual = FindLongestRepeat(Input[0],Input[1])
        print (len(Expected[0]),len(Actual))
        print (Expected[0])
        print (Actual)        
 
    if args.rosalind:
        pass
    
    elapsed = time.time()-start
    minutes = int(elapsed/60)
    seconds = elapsed-60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')