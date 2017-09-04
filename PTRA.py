from Bio.Seq import translate
import sys

with open ('c:/Users/Weka/Downloads/rosalind_ptra(2).txt') as f:
    i=0
    for line in f:
        print (line.strip())
        if i==0:
            coding_dna=line.strip()
        if i==1:
            out=line.strip()
            for table in [1,2,3,4,5,6,9,10,11,12,13,14,15]:
                print (table)
                translated=translate(coding_dna,table=table,to_stop=True)
#                if len(translated)==len(out):
                if translated==out:
                    print (table,translated)
                    sys.exit()
                print ('Not matched: '+ translated)
        i+=1
       