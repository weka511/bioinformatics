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

import re
from urllib.request import urlopen
from urllib.parse import urljoin
def mprt(url = 'http://www.uniprot.org/uniprot/',ID='B5ZC00',motif_pattern='N[A-O,Q-Z][ST][A-O,Q-Z]'):
      resource = urlopen(urljoin(url,ID))
      content  = resource.read().decode(resource.headers.get_content_charset())
      re_motif = re.compile(motif_pattern)
      motifs   = {}
      pos      = content.find('Sequence status')
      while True:
            matched = re_motif.search(content,pos=pos)
            if matched:
                  motif         = matched.string[matched.start(0):matched.end(0)]
                  pos           = matched.end(0)
                  while content[pos].isalpha():
                        pos+=1
                  start = matched.start(0)
                  while content[start-1].isalpha():
                        start-=1                  
                  L = 59

                  try: 
                        length = matched.end(0) - matched.start(0)
                        n = int(content[start-L:pos-L])
                        offset = pos - matched.end(0)
                        p = n -length -offset +1
                        #print (f'{motif}--{n}--{content[start-L:pos-L]}--{p}')
                        motifs[motif] = (matched.start(0),start,matched.end(0),pos,p)
                  except:
                        ValueError
            else:
                  break
      return motifs


      
if __name__=='__main__':
      for ID in ['A2Z669','B5ZC00','P07204_TRBM_HUMAN','P20840_SAG1_YEAST']:
            motifs = mprt(ID=ID)
            if len(motifs)>0:
                  print (ID)
                  positions = [pos for _,_,_,_,pos in motifs.values()]
                  positions.sort()
                  print (' '.join([str(p) for p in positions]))