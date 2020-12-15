# Copyright (C) 2017-2020 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

# This file contains a collection of functions to parse FASTA format data

# FastaContent
#
# This class is used to parse FASTA format data from a text string

class FastaContent(object):
    # parse text and build data structures
    def __init__(self,file_text):
        self.pairs = []             # Description + data
        nn         = ''             # Description
        text       = ''             # data
        for line in file_text:
            line = line.strip()
            if line.startswith('>'):     # Description
                if len(nn)>0:
                    self.pairs.append((nn,text))
                nn   = line[1:]
                text = ''
            else:                         # Data
                text=text+line
                
        if len(nn)>0:
            self.pairs.append((nn,text))
    
    # Retrieve items (pairs) by index       
    def __getitem__(self,index):
        return self.pairs[index]
    
    # Number of entries
    def __len__(self):
        return len(self.pairs)
    
    # Used if caller wants a simple list seq,value, seq,value,...
    
    def to_list(self):
        List = []
        for a,b in self.pairs:
            List.append(a)
            List.append(b)
        return List

# FastaFile
#
# This class is used to parse Fasta data from a file
class FastaFile(FastaContent):
    # parse file and build data structures
    def __init__(self,name):
        with open(name,'r') as file_text:
            FastaContent.__init__(self,file_text)
  
# fasta_out
#
# Generator, used to output one key value pair to a file in FASTA format
#
# Parameters:
#      key           Omit '>', as this will be added
#      value         Value string
#      max_length    Used to split lines so none exceeds this length

def fasta_out(key,value,max_length=80):
    yield f'>{key}'
    remainder = value
    while len(remainder)>0:
        yield remainder[0:max_length]
        remainder = remainder[max_length:]
            