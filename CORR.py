'''
 CORR Error Correction in Reads
 
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



def corr(reads):
	def create_read_match_counts(reads):
		read_match_counts={}
	
		for key,read in reads:
			if read in read_match_counts:
				read_match_counts[read]+=1
			else:
				read_match_counts[read]=1
		return read_match_counts

	def create_matches(read_match_counts):
		matches={}
		for key,count in read_match_counts.items():
			if count>=2:
				matches[key]=count
		return matches

	def create_unmatched(read_match_counts,matches):
		unmatched=[]		
		for key,count in read_match_counts.items():
			if count==1 and not key in matches:
				r.revc_key=r.revc(key)
				if r.revc_key in read_match_counts:
					matches[key]=read_match_counts[r.revc_key]+1
					matches[r.revc_key]=read_match_counts[r.revc_key]+1
				else:
					unmatched.append(key)
		return unmatched
	
	def find_neighbours(key,matches):
		neighbours=[]
		for i in range(len(key)):
			for letter in ['A','C','G','T']:
				if letter!= key[i]:
					key_1=key[:i] + letter+key[i+1:]
					if key_1 in matches:
						neighbours.append(key_1)
		return neighbours
	
	def fixup(unmatched,matches):
		pairs=[]
		for read in unmatched:
			neighbours=find_neighbours(read,matches)
			r.revc_neighbours= find_neighbours(r.revc(read),matches)
			if len(neighbours)==1 and len(r.revc_neighbours)==0:
				pairs.append((read,neighbours[0]))
			if len(neighbours)==0 and len(r.revc_neighbours)==1:
				pairs.append((read,r.revc(neighbours[0])))
		return pairs
	
	read_match_counts = create_read_match_counts(reads)
	matches=create_matches(read_match_counts)
	unmatched=create_unmatched(read_match_counts,matches)
	return fixup(unmatched,matches)



if __name__=='__main__':    
	from Bio import SeqIO
	def combine(seq_record):
		return (seq_record.id,str(seq_record.seq))        
	for a,b in corr([combine(seq_record) for seq_record in SeqIO.parse("c:/Users/Weka/Downloads/rosalind_corr(6).txt", "fasta")]):
		print ('{0}->{1}'.format(a,b))