# Copyright (C) 2019 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# BA9C Construct the Suffix Tree of a String 

# ConstructModifiedSuffixTrie
#
# Construct Modified Suffix Trie, as described in 
# Charging Station Bioinformatics Algorithms Vol II page 165
#
# See also https://sandipanweb.wordpress.com/2017/05/10/suffix-tree-construction-and-the-longest-repeated-substring-problem-in-python/

    


def ConstructSuffixTreeEdges(s):
    def explore(suffixes):
        if len(suffixes)==0: return
        prefixes = list(set([s[0] for s in suffixes]))
        if len(prefixes)==len(suffixes):
            for s in suffixes:
                Edges.append(s)
        else:
            if len(prefixes)==1:
                for k in range(2,min([len(s) for s in suffixes])):
                    extended_prefixes = list(set([s[0:k] for s in suffixes]))
                    if len(extended_prefixes)==1:
                        prefixes = extended_prefixes
                    else:
                        break
                    
            for p in prefixes:
                subset = [s[len(prefixes[0]):] for s in suffixes if s[0:len(prefixes[0])]==p ]
                if len(subset)>1:
                    Edges.append(p)
                    explore(subset)
                elif len(subset)==1:
                    Edges.append(p+subset[0])
                #else:
                    #Edges.append(p)
                    
     
    next_node = 0                
    Edges     = []
    explore (sorted([s[i:] for i in range(len(s))]))
    return Edges

if __name__=='__main__':
    #print('---')
    #for edge in ConstructSuffixTreeEdges('ATAAATG$'):
        #print(edge)
    #print('---')
    Edges = ConstructSuffixTreeEdges('ATCTACCAGCAGTGAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$')   
  
    Expected = [
        '$',
        'ATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'GCCGCCGCTGG$',
        'C',
        'AAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'A',
        'GG',
        'GACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'TGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'A',
        'CAGTGAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'TAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'AG',
        'CTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'C',
        'GCCGCCGCTGG$',
        'TCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'C',
        'ACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CAGTGAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'AGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'A',
        'CTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'GTGTATAACGCCGCCGCTGG$',
        'G',
        'AAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'GAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'G',
        'T',
        'AACGCCGCCGCTGG$',
        'CTACCAGCAGTGAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'GGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'TGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'A',
        'T',
        'G',
        'ACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CAGTGAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'AAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'GAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'G',
        'T',
        'TGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'A',
        'CAGTGAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'TAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'AG',
        'CTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CGCTGG$',
        'TGG$',
        'GC',
        'TCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'C',
        'ATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CGCTGG$',
        'TGG$',
        'CGC',
        'TGG$',
        'C',
        'AGGGTGTATAACGCCGCCGCTGG$',
        'TCGTAGGGTGTATAACGCCGCCGCTGG$',
        'G',
        'T',
        'ACCAGCAGTGAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'ATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'TTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CG',
        'GG$',
        'TACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'C',
        'T',
        '$',
        'CATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'GGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'A',
        'CAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'TCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'C',
        'GGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'TGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'A',
        'AGTGAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CGCTGG$',
        'TGG$',
        'CGC',
        'GG$',
        'TACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'C',
        'T',
        '$',
        'AGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'GGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'A',
        'CTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'AGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'TGTATAACGCCGCCGCTGG$',
        'G',
        'TGTATAACGCCGCCGCTGG$',
        'G',
        'AGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'GGGTGTATAACGCCGCCGCTGG$',
        'TAACGCCGCCGCTGG$',
        'A',
        'AACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'ATAACGCCGCCGCTGG$',
        'TACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'G',
        'T',
        'ACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CGTAGGGTGTATAACGCCGCCGCTGG$',
        'G',
        'T',
        'T',
        'CGCCGCCGCTGG$',
        'GGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'A',
        'AGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'AGCAGTGAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'C',
        'C',
        'GGGTGTATAACGCCGCCGCTGG$',
        'TAACGCCGCCGCTGG$',
        'A',
        'ATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'AGGGTGTATAACGCCGCCGCTGG$',
        'TCGTAGGGTGTATAACGCCGCCGCTGG$',
        'G',
        'T',
        'TACCAGCAGTGAACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'C',
        'AACATGGGAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        '$',
        'GAGGACCAGTAAGGAAGGCTTACCCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'G',
        'ATAACGCCGCCGCTGG$',
        'GTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'TACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'G',
        'T',
        'AGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'CCTCGATGTGTTACAGACTCGTTCGTAGGGTGTATAACGCCGCCGCTGG$',
        'AC',
        'T',
        'T',
        'CGTAGGGTGTATAACGCCGCCGCTGG$'
    ]
    
    print ('Expected = {0}, actual={1}'.format(len(Expected),len(Edges)))
    expected = iter(sorted(Expected))
    edges    = iter(sorted(Edges))
    exp = next(expected)
    ed = next(edges)
    while exp != '-' and ed !='-':
        if exp<ed:
            print('{0},{1}'.format(exp,'-'))
            exp = next(expected,'-')           
        elif ed<exp:
            print('{0},{1}'.format('-',ed))
            ed = next(edges,'-') 
        else:
            exp = next(expected,'-')
            ed = next(edges,'-')            
            