from Bio import Entrez,SeqIO
inputs = 'JN698988 NM_001271262 JN573266 JQ342169 NM_001133698 JX462669 NM_204821 BT149867 JX472277'

Entrez.email = "simon@greenweaves.nz"
handle = Entrez.efetch(db="nucleotide", id=[inputs.replace(' ',', ')], rettype="fasta")
records = list (SeqIO.parse(handle, "fasta"))
best=None
first=True
for record in records:
    if first or len(record)<len(best):
        best = record
        first = False
print(best.format("fasta"))

    