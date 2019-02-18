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
# EDIT 	Edit Distance http://rosalind.info/problems/edit/

def edit(s,t,indel_cost=1,replace_cost=lambda a,b: 1,show_matrix=False):
    
    def dynamic_programming(s,t):
        matrix=[[0 for j in range(len(t)+1)] for i in range(len(s)+1)]
    
        for j in range(len(t)+1):
            matrix[0][j]=j
        for i in range(len(s)+1):
            matrix[i][0]=i

        for i in range(1,len(s)+1):
            for j in range(1,len(t)+1):
                matrix[i][j] = min(
                    matrix[i-1][j]   + indel_cost,
                    matrix[i][j-1]   + indel_cost,
                    matrix[i-1][j-1] + (0 if s[i-1]==t[j-1] else replace_cost(s[i-1],t[j-1])))
                    
        if show_matrix:
            for i in range(0,len(s)+1):
                ii = len(matrix)-i-1
                print (s[ii] if i>0 else '#',matrix[ii])
            print (' ',['#']+t)
        
        return matrix[len(s)][len(t)],matrix
            
    return dynamic_programming([s0 for s0 in s], [t0 for t0 in t])

if __name__=='__main__':
    #from Bio import SeqIO
    #inFile = open(r'C:\Users\Simon\Downloads\rosalind_ba5g(1).txt','r')
    #strings = []
    #for record in SeqIO.parse(inFile,'fasta'):
        #print (record.id)
        #print (str(record.seq))
    #strings.append(str(record.seq))
    #d,_=edit(strings[0],strings[1])
    #print(d)
    with open(r'C:\Users\Simon\Downloads\rosalind_ba5g(2).txt','r') as f:
        strings = []
        for line in f:
            strings.append(line.strip())
            print (line.strip())
        d,_=edit(strings[0],strings[1])
        print(d)
            
        