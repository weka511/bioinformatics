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
#   MPRT

import re
from urllib.request import urlopen
from urllib.parse import urljoin
def mprt(url = 'http://www.uniprot.org/uniprot/',ID='B5ZC00',motif_pattenr='N[A-O,Q-Z][ST][A-O,Q-Z]'):
      resource = urlopen(urljoin(url,ID))
      content  = resource.read().decode(resource.headers.get_content_charset())
      re_motif = re.compile(motif_pattenr)
      motifs   = set()
      pos = 0
      while True:
            matched = re_motif.search(content,pos=pos)
            if matched:
                  motif = matched.string[matched.start(0):matched.end(0)]
                  motifs.add(motif)
                  pos = matched.end(0)
            else:
                  break
      print (motifs)
      
if __name__=='__main__':
      mprt()