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

# LCSQ 	Finding a Shared Spliced Motif 

# See http://www.cs.cmu.edu/afs/cs/academic/class/15451-s15/LectureNotes/lecture04.pdf, and
# https://www.ics.uci.edu/~goodrich/teach/cs260P/notes/LCS.pdf

def dynamic(s,t,extract):
    matrix=[[0 for j in t] for i in s]
    def get_matrix(i,j):
        return matrix[i][j] if i>=0 and j>=0 else 0            
    for i in range(len(s)):
        for j in range(len(t)):
            matrix[i][j] = (get_matrix(i-1,j-1) + 1) if s[i]==t[j] else max(get_matrix(i-1,j), get_matrix(i,j-1))
    return extract(s,t,get_matrix)
      
def lcsq(s,t):
    def extract(s,t,get_matrix):        
        i = len(s)-1
        j = len(t)-1
        seq = []
        while True:
            if get_matrix(i,j-1) == get_matrix(i,j): j -= 1
            elif get_matrix(i-1,j) == get_matrix(i,j): i -= 1
            else:
                assert s[i]==t[j]
                seq.append(s[i])
                i -= 1
                j -= 1
            if  get_matrix(i,j)==0:
                return seq[::-1]    
    return ''.join(dynamic([s0 for s0 in s], [t0 for t0 in t],extract))

if __name__=='__main__':
    print (lcsq('GATTCAGAGAATGTTTTCGTAAGCTGTAATTCGGTTTGTATTTGCAACGAACGATTGACA'
                'GAAAGAGGGATATAACGGCAGACGTACGGAAAGAACTCGAGCCATGCCATCATACGAACT'
                'TCGCGTAAGTATTGGGGCCTTAGGGTGTCATGGCCCAGCTATGTACAAAATGATCATGTC'
                'GTTGACAACCTTATGAGAGTAGATGACCTCTAACCACCCGAAAGGGCTACGGTTAGAACT'
                'AGCGAAATATATTAGCCTCTTACCTGCTGGGGATAGGCATAGTTCACAACGTCTCCGGCC'
                'GCGATTTGATTACCTTACGCTGGGAGGCTTCCCCCCTCGGGATTCGGGGAGACTCATTAC'
                'AACGGTGCACAGAGGAATGATTTTCCTAAGATATGTACTCTAACCAGCAAGATGGGCATA'
                'CATAGCCGGCTTTGACTGTAATTATGTGTGTTACATGAGCCGTATCAGCAGTTGGTTGCG'
                'ATGCCAGAACGAGTCAGTATGATAACTCCCGGCATAACGGAACCCAAACCGGTCGCTACA'
                'TCTCGTGTATAAAGGGCGTACCGCCCCCTCCCCAGACAGGCCGCACCATAGGGTCCTTAG'
                'TGTGATCATGGCGGAGCACCGGATTGGTGGGTTGGCATATCTCCGCGCGTGCGATTTCTG'
                'TCCCTCATTCGACTCCAACTACCCCGACTTGCGTCGTAAAGGTAAACATTCCGGTCTCAT'
                'ATGGAGGCTGAGGAGCTTTAGAGAAGAAAAAAGCCCTATAACGTTATCAGGGTCCTCTTT'
                'CGTACGAGTTCCAAAATTTTCCTTCCGCCGGCTTGTTCAAATGATGTATGGAGATTCTTC'
                'CCTATAGGCGGACGTTGAGTATAATGCCCGATGTAGAAGTCCAGGGTGGGTATCATAATA'
                'GTATGCCCGGTGACAAATAACGAATGTCTATGCGAGTCACGGCAGTAA',
                'ATGCTATATTGTGACATACTCACACGTTTACGCTATGTCCTGAAAATTGTGCGGGTTGCC'
                'AACTCGCCGGGCTTAGACAGCGTAACTTGCGCTATCCCATCGTTGTTCGATGAACTACCT'
                'CACTAAAGTCATTTATAGTTTAGCTCGCCACTGATCTTCGGTCGACGGTTGATCACACGT'
                'CGAGGCATTTTCCCACACGAGGGATCCGGATCAACTATCTCCGGTAGATAGAGTAAGTAG'
                'GCCCTTTAGCTCGACACACCTTGGTACGATTTGCCGCAATCTGCTCGAGTCTTTTTTTTC'
                'TTCTGCAGAACGGTTGGAAGCATTTTCGAGTTCCTCTAGGTCCCCGACTACAGATTTCAG'
                'CGATCTGATAGGTAACGGTCCAGATCAGACCTCGAACGCTCTAACAATATACAGAGTTTA'
                'CTCTAAACCAGCCGTCTGAGTGATTATTTGGAGAACCTACTGCGACCGTGCAAATGCAGG'
                'CGATACTGGTGGTCCATGGAAGGACGCGTCCCTCGGCGAAATACATATTTTCAGCCCAAC'
                'CGCCGCGTAGGATGAGACTTTTAGAGTCGCACTAACAGTCTCTAGGGAATAGGGAGTCAA'
                'GTGCCAGGCCAGCTCGTTGACTTGGTGCGGGACCCTGACACCAAAGTTGTGTACTAGTCA'
                'AATCACGGGGCGCAGGAAACCTCACCCTTTATTACGGTACATTGAACTCATCGTAATTAC'
                'CATGGGATTAACAGGGCGCACCATTGGCTTCTTCCTAGCCCTCACCCTGGTTCGTGAGAA'
                'TCCGGGTGAGCTTTCTTGTTACACCAGTCACAGATAACAACATAAGTGTTTGGGGTTCCC'
                'GGAATTC'))
    