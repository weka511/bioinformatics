def revc(dna):
	return dna.translate({
	    ord('A'): 'T',
	    ord('C'): 'G',
	    ord('G'):'C', 
	    ord('T'): 'A'})[::-1]

def corr(strings):
	def fixup1(key):
		for i in range(len(key)):
			for letter in ['A','C','G','T']:
				if letter!= key[i]:
					key_1=key[:i] + letter+key[i+1:]
					if key_1 in dnas:
						matches.append(key_1)
						return key_1
		return None
	def fixup(key,revc_key):
		fixit=fixup1(key)
		if fixit==None:
			revc_fixed=fixup1(revc_key)
			return key,revc(revc_fixed)
		else:
			return key,fixit
	dnas={}
	pairs=[]
	matches=[]
	for key,string in strings:
		if string in dnas:
			dnas[string]+=1
		else:
			dnas[string]=1
	for key,count in dnas.items():
		if count==1 and not key in matches:
			revc_key=revc(key)
			if not revc_key in dnas:
				pairs.append(fixup(key,revc_key))
	return pairs

if __name__=='__main__':    
	from Bio import SeqIO
	def combine(seq_record):
		return (seq_record.id,str(seq_record.seq))        
	#print (long([combine(seq_record) for seq_record in SeqIO.parse("c:/Users/Weka/Downloads/rosalind_long.txt", "fasta")]))
	for a,b in corr([combine(seq_record) for seq_record in SeqIO.parse("c:/Users/Weka/Downloads/rosalind_corr(4).txt", "fasta")]):
	#for a,b in corr([combine(seq_record) for seq_record in SeqIO.parse("CORR.txt", "fasta")]):
		print ('{0}->{1}'.format(a,b))