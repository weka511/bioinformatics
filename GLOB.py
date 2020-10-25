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
#    GLOB Global Alignment with Scoring Matrix

if __name__=='__main__':
    from  align import highest_scoring_global_alignment
    from Bio import SeqIO
    inFile = open(r'C:\Users\Simon\Downloads\rosalind_glob.txt','r')
    strings = []
    for record in SeqIO.parse(inFile,'fasta'):
        print (record.id)
        print (str(record.seq))
        strings.append(str(record.seq))
    d,s,t=highest_scoring_global_alignment(strings[0],strings[1])
    print(d)
    print (s)
    print (t)    
    print (highest_scoring_global_alignment('PLEASANTLY','MEANLY'))