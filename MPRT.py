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

def remove_white_space(motif):
      pos_space = motif.find(' ')
      return ''.join(motif[i] for i in range(len(motif)) if motif[i].isalpha() and motif[i].isupper()),pos_space

def mprt(url = 'http://www.uniprot.org/uniprot/',ID='B5ZC00',motif_pattern='N[\s]*[^P][\s]*[ST][\s]*[^P]'):
      L        = 59
      resource = urlopen(urljoin(url,ID))
      content  = resource.read().decode(resource.headers.get_content_charset())
      re_motif = re.compile(motif_pattern)
      motifs   = {}
      pos      = content.find('Sequence status')
      final    = content.find('Sequence databases')
      while True:
            matched = re_motif.search(content,pos=pos)
            if matched:
                  motif,pos_space = remove_white_space(matched.string[matched.start(0):matched.end(0)])
                  if matched.start(0)>final:
                        break
                  pos=matched.end(0)
                  if len(motif)!=4:
                        continue
                  #print (motif, matched.string[matched.start(0):matched.end(0)])
                  seq_loc_start = matched.start(0) - L
                  seq_loc_end   = matched.start(0) - L
                  while content[seq_loc_start].isdigit():
                        seq_loc_start-=1
                  if seq_loc_start!=' ':
                        seq_loc_start+=1  
                  while not content[seq_loc_end].isdigit():
                        seq_loc_end+=1
                  while content[seq_loc_end].isdigit():
                        seq_loc_end+=1
                  try:
                        xx=content[seq_loc_start:seq_loc_end]
                        loc_upper_bound = int(content[seq_loc_start:seq_loc_end].strip())
                        seq = loc_upper_bound - ( seq_loc_end-(matched.end(0)-L)) -len(motif) +1
                        if pos_space>-1:
                              seq-=(len(motif)-pos_space)
                        #print(motif,loc_upper_bound, seq_loc_start - matched.start(0)+L, 
                              #matched.end(0)-L, pos_space, seq_loc_end, seq_loc_end-(matched.end(0)-L),seq)
         
                        motifs[motif] = seq
                  except ValueError:
                        break
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