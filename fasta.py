#!/usr/bin/env python
# Copyright (C) 2017-2024 Greenweaves Software Limited

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

''' This file contains a collection of functions to parse FASTA format data'''


class FastaContent(object):
    '''
    FastaContent

    This class is used to parse FASTA format data from a list of text strings,
    each representing one line from a file
    '''
    def __init__(self,file_text):
        '''parse text and build data structures'''
        self.pairs = []             # Description + data -- all text so fae
        nn = ''             # Description        -- current
        text = ''             # data               -- current
        for line in file_text:
            line = line.strip()
            if line.startswith('>'):     # Description
                if len(nn)>0:
                    self.pairs.append((nn,text))
                nn = line[1:]
                text = ''
            else:                         # Data
                text = text + line

        if len(nn)>0:
            self.pairs.append((nn,text))


    def __getitem__(self,index):
        '''Retrieve items (pairs) by index'''
        return self.pairs[index]


    def __len__(self):
        '''Number of entries'''
        return len(self.pairs)

    def to_list(self):
        '''
        to_list

        Used if caller wants a simple list seq,value, seq,value,...
        '''
        List = []
        for a,b in self.pairs:
            List.append(a)
            List.append(b)
        return List


    def to_dict(self):
        '''
        Used if caller wants a dictionary
        '''
        Dict = {}
        for a,b in self.pairs:
            Dict[a]=b
        return Dict


class FastaFile(FastaContent):
    '''
    FastaFile

    This class is used to parse Fasta data from a file
    '''
    def __init__(self,name):
        ''' parse file and build data structures'''
        with open(name,'r') as file_text:
            FastaContent.__init__(self,file_text)


def fasta_out(key,value,
              max_length=800000000000):
    '''
    fasta_out

    Generator, used to output one key value pair to a file in FASTA format

    Parameters:
        key           Omit '>', as this will be added
        value         Value string
        max_length    Used to split lines so none exceeds this length
    '''
    yield f'>{key}'
    remainder = value
    while len(remainder)>0:
        yield remainder[0:max_length]
        remainder = remainder[max_length:]
