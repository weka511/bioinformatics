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

import os
import re
from urllib.request import urlopen
from urllib.parse import urljoin
import reference_tables as rt

# read_uniprot_as_fasta
#
# Read data from uniprot site as fasta file     

def read_uniprot_as_fasta(url = 'http://www.uniprot.org/uniprot/',ID='B5ZC00'):
      resource = urlopen(urljoin(url,ID+'.FASTA'))
      return resource.read().decode("utf-8", "ignore")

# get_protein_sequence
#
# Extract the portion of the fasta stringthat represents protein sequence

def get_protein_sequence(fasta):
      pos         = 0
      re_protein  = rt.get_re_protein(min_length=8)
      matched     = re_protein.search(fasta,pos=pos)
      amino_acids = []
      while matched:
            amino_acids.append(matched.string[matched.start(0):matched.end(0)])
            pos        = matched.end(0)+1
            re_protein = rt.get_re_protein(min_length=1)
            matched    = re_protein.search(fasta,pos=pos)
      return ''.join(amino_acids) 

# mprt
#
# Input: At most 15 UniProt Protein Database access IDs.
#
# Return: For each protein possessing the N-glycosylation motif, output its given access ID followed
#         by a list of locations in the protein string where the motif can be found.

def mprt(url = 'http://www.uniprot.org/uniprot/',ID='B5ZC00',motif_pattern='N[\s]*[^P][\s]*[ST][\s]*[^P]'):
      fasta            = read_uniprot_as_fasta(url=url,ID=ID)
      protein_sequence = get_protein_sequence(fasta)
      start_protein    = fasta.find('\n')
      re_motif         = re.compile(motif_pattern)
      motifs           = {}
      matched          = re_motif.search(protein_sequence,pos=0)
      while matched:
            motif         = matched.string[matched.start(0):matched.end(0)]
            pos           = matched.start(0) + 1
            motifs[motif] = pos
            matched       = re_motif.search(protein_sequence,pos=pos)

      return motifs

def read_data(file,path=r'C:\Users\Simon\bioinformatics\data'):
      with open (os.path.join(path,file)) as f:
            return [line.strip() for line in f]
      
if __name__=='__main__':
      import argparse
      parser = argparse.ArgumentParser('MPRT')
      parser.add_argument('file', nargs='?')
      parser.add_argument('--path', default = r'C:\Users\Simon\bioinformatics\data')
      parser.add_argument('--output')
      args=parser.parse_args()
      IDs = read_data(args.file,path=args.path) if args.file != None else [
            'A2Z669',
            'B5ZC00',
            'P07204_TRBM_HUMAN',
            'P20840_SAG1_YEAST'
            ]
      if args.output:
            with open (args.output,'w') as out:
                  for ID in IDs:    
                        motifs = mprt(ID=ID)
                        if len(motifs)>0:
                              out.write (f'{ID}\n')
                              locations = [seq for seq in motifs.values()]
                              locations.sort()
                              line = ' '.join([str(p) for p in locations])
                              out.write (f'{line}\n')                  
      else:            
            for ID in IDs:    
                  motifs = mprt(ID=ID)
                  if len(motifs)>0:
                        print (ID)
                        locations = [seq for seq in motifs.values()]
                        locations.sort()
                        print (' '.join([str(p) for p in locations]))