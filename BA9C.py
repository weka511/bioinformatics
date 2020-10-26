#  Copyright (C) 2019 Greenweaves Software Limited
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
# BA9C Construct the Suffix Tree of a String 

import argparse
from snp import ConstructSuffixTreeEdges    

def check(Edges,Expected):
    mismatches = 0
    for a,b in zip(sorted(Edges),sorted(Expected)):
        if a!=b:
            mismatches+=1
            print (f'Expected {b}, was {a}')
    print(f'{0} mismatches')
    
if __name__=='__main__':
    parser = argparse.ArgumentParser('BA9C Construct the Suffix Tree of a String ')
    parser.add_argument('--sample', default=False, action='store_true', help='process sample dataset')
    args = parser.parse_args()
    if args.sample:
        Edges = ConstructSuffixTreeEdges('ATAAATG$')
        print('---')
        for edge in Edges:
            print(edge)
        print('---')
        Expected = [
            'AAATG$',
            'G$',
            'T',
            'ATG$',
            'TG$',
            'A',
            'A',
            'AAATG$',
            'G$',
            'T',
            'G$',
            '$'        
        ]
        check(Edges,Expected)
    else:
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
        exp      = next(expected)
        ed       = next(edges)
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
            