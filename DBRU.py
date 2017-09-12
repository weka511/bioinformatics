'''
 DBRU Constructing a De Bruijn Graph

 Copyright (C) 2017 Greenweaves Software Pty Ltd

 This is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This software is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

 As is the case with point mutations, the most common type of sequencing
 error occurs when a single nucleotide from a read is interpreted incorrectly.

 Given: A collection of up to 1000 reads of equal length (at most 50 bp) in
 FASTA format. Some of these reads were generated with a single-nucleotide error.
 For each read s in the dataset, one of the following applies:
    s was correctly sequenced and appears in the dataset at least twice
	   (possibly as a reverse complement);
    s is incorrect, it appears in the dataset exactly once, and its
	  Hamming distance is 1 with respect to exactly one correct read 
	  in the dataset (or its reverse complement).

 Return: A list of all corrections in the form "[old read]->[new read]". 
 (Each correction must be a single symbol substitution, and you may return the corrections in any order.)
'''
import rosalind as r

def dbru(S,include_revc=True):
    def union(S):
        U=set(S)
        for s in S:
            revc=r.revc(s)
            if include_revc and not revc in U:
                U.add(revc)
        return U
 
    return [(e[0:-1],e[1:]) for e in union(S)]

if __name__=='__main__':
    A=[]
    with open('c:/Users/Weka/Downloads/rosalind_dbru(1).txt') as f:
        for line in f:
            A.append(line.strip())       
    #for a,b in dbru(['TGAT','CATG','TCAT','ATGC','CATC','CATC']):
    for a,b in dbru(A):
        print('({0}, {1})'.format(a,b))
