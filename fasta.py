# Copyright (C) 2020 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>

# This file contains a collection of functions to solve the problems
# at rosalind.info.

# FastaContent
#
# This class is used to parse FASTA format data from a text string

class FastaContent(object):
    def __init__(self,file_text):
        self.pairs = []
        nn         = ''
        text       = ''    
        for line in file_text:
            line = line.strip()
            if line.startswith('>'):
                if len(nn)>0:
                    self.pairs.append((nn,text))
                nn   = line[1:]
                text = ''
            else:
                text=text+line
                
        if len(nn)>0:
            self.pairs.append((nn,text))
            
    def __getitem__(self,index):
        return self.pairs[index]
    
    def __len__(self):
        return len(self.pairs)

# FastaFile
#
# This class is used to parse Fasta data from a file
class FastaFile(FastaContent):
    def __init__(self,name):
        with open(name,'r') as file_text:
            FastaContent.__init__(self,file_text)
            