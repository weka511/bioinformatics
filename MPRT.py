#    Copyright (C) 2020 Simon Crase
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
#   MPRT Finding a Protein Motif
# 
#  N{P}[ST]{P}.
import re
from urllib.request import urlopen
from urllib.parse import urljoin
import reference_tables as rt

      
def read_unipro_as_fasta(url = 'http://www.uniprot.org/uniprot/',ID='B5ZC00'):
      resource = urlopen(urljoin(url,ID+'.FASTA'))
      return resource.read().decode("utf-8", "ignore")

def get_protein_sequence(fasta):
      pos=0
      re_protein = rt.get_re_protein(min_length=8)
      matched = re_protein.search(fasta,pos=pos)
      amino_acids=[]
      while matched:
            amino_acids.append(matched.string[matched.start(0):matched.end(0)])
            #print (pos, matched.string[matched.start(0):matched.end(0)])
            pos = matched.end(0)+1
            re_protein = rt.get_re_protein(min_length=1)
            matched = re_protein.search(fasta,pos=pos)
      return ''.join(amino_acids)            
def mprt(url = 'http://www.uniprot.org/uniprot/',ID='B5ZC00',motif_pattern='N[\s]*[^P][\s]*[ST][\s]*[^P]'):
      fasta = read_unipro_as_fasta(url=url,ID=ID)
      protein_sequence = get_protein_sequence(fasta)
      start_protein = fasta.find('\n')
      re_motif = re.compile(motif_pattern)
      motifs   = {}
      pos = 0
      while True:
            matched = re_motif.search(protein_sequence,pos=pos)
            if matched:
                  motif = matched.string[matched.start(0):matched.end(0)]
                  pos = matched.start(0) + 1
                  matched = re_motif.search(protein_sequence,pos=pos)
                  motifs[motif] = pos
            else:
                  break
      return motifs


      
if __name__=='__main__':
      for ID in [
            'A2Z669',
            'B5ZC00',
            'P07204_TRBM_HUMAN',
            'P20840_SAG1_YEAST'
            ]:
            #'P04141_CSF2_HUMAN',
            #'P01190_COLI_BOVIN',
            #'P07204_TRBM_HUMAN',
            #'Q8WW18',
            #'Q3B391',
            #'P0AF66',
            #'Q5FTZ8',
            #'P80069_A45K_MYCBO',
            #'P37803',
            #'Q14ID0',
            #'P80195_MPP3_BOVIN',
            #'P11279_LMP1_HUMAN',
            #'A5GIU0',
            #'A7Z201'               
            motifs = mprt(ID=ID)
            if len(motifs)>0:
                  print (ID)
                  locations = [seq for seq in motifs.values()]
                  locations.sort()
                  print (' '.join([str(p) for p in locations]))