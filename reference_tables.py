# Copyright (C) 2015-2020 Greenweaves Software Limited

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

import re
bases='ACGT'

codon_table={
     'UUU': 'F',      'CUU': 'L',      'AUU': 'I',      'GUU': 'V',
     'UUC': 'F',      'CUC': 'L',      'AUC': 'I',      'GUC': 'V',
     'UUA': 'L',      'CUA': 'L',      'AUA': 'I',      'GUA': 'V',
     'UUG': 'L',      'CUG': 'L',      'AUG': 'M',      'GUG': 'V',
     'UCU': 'S',      'CCU': 'P',      'ACU': 'T',      'GCU': 'A',
     'UCC': 'S',      'CCC': 'P',      'ACC': 'T',      'GCC': 'A',
     'UCA': 'S',      'CCA': 'P',      'ACA': 'T',      'GCA': 'A',
     'UCG': 'S',      'CCG': 'P',      'ACG': 'T',      'GCG': 'A',
     'UAU': 'Y',      'CAU': 'H',      'AAU': 'N',      'GAU': 'D',
     'UAC': 'Y',      'CAC': 'H',      'AAC': 'N',      'GAC': 'D',
     'UAA': ';',      'CAA': 'Q',      'AAA': 'K',      'GAA': 'E',
     'UAG': ';',      'CAG': 'Q',      'AAG': 'K',      'GAG': 'E',
     'UGU': 'C',      'CGU': 'R',      'AGU': 'S',      'GGU': 'G',
     'UGC': 'C',      'CGC': 'R',      'AGC': 'S',      'GGC': 'G',
     'UGA': ';',      'CGA': 'R',      'AGA': 'R',      'GGA': 'G',
     'UGG': 'W',      'CGG': 'R',      'AGG': 'R',      'GGG': 'G'
 }

integer_masses={
    'G': 57,
    'A': 71,
    'S': 87,
    'P': 97,
    'V': 99,
    'T': 101,
    'C': 103,
    'I': 113,
    'L': 113,
    'N': 114,
    'D': 115,
    'K': 128,
    'Q': 128,
    'E': 129,
    'M': 131,
    'H': 137,
    'F': 147,
    'R': 156,
    'Y': 163,
    'W': 186
}

# AminoAcid
#
# This class represents the aspects of amino acids that are relevant to 
# Mass Spectroscopy

#   Attributes:
#     name         Name of amino acid, e.g. Alanine
#     short        Single letter code, e.g. A for Alanine
#     abbrev       Three digit code, e.g. Ala
#     mon_mass     Monoisotpic mass (Daltons)
#     average_mass Average mass (Daltons)

class AminoAcid:
    
    def __init__(self,name,short,abbrev,mon_mass,average_mass):
        self.name         = name
        self.short        = short
        self.abbrev       = abbrev
        self.mon_mass     = mon_mass
        self.average_mass = average_mass
        
    def __str__(self):
        return '%(name)s %(short)s %(abbrev)s %(int_mass)d %(mon_mass)f %(average_mass)f'%\
               {
                   'name'         : self.name,
                   'short'        : self.short,
                   'abbrev'       : self.abbrev,
                   'mon_mass'     : self.mon_mass,
                   'average_mass' : self.average_mass,
                   'int_mass'     : self.asInteger()
               }
    
    def asInteger(self):
        return int(self.mon_mass)
 
# amino_acids
#
# Lookup table for amino acids, from
# https://en.wikipedia.org/wiki/Proteinogenic_amino_acid#Mass_spectrometry

amino_acids={
    'A': AminoAcid('Alanine',        'A', 'Ala',  71.03711,  71.0788),
    'C': AminoAcid('Cysteine',       'C', 'Cys', 103.00919, 103.1388),
    'D': AminoAcid('Aspartic acid',  'D', 'Asp', 115.02694, 115.0886),
    'E': AminoAcid('Glutamic acid',  'E', 'Glu', 129.04259, 129.1155),
    'F': AminoAcid('Phenylalanine',  'F', 'Phe', 147.06841, 147.1766),
    'G': AminoAcid('Glycine',        'G', 'Gly',  57.02146,  57.0519),
    'H': AminoAcid('Histidine',      'H', 'His', 137.05891, 137.1411),
    'I': AminoAcid('Isoleucine',     'I', 'Ile', 113.08406, 113.1594),
    'K': AminoAcid('Lysine',         'K', 'Lys', 128.09496, 128.1741),
    'L': AminoAcid('Leucine',        'L', 'Leu', 113.08406, 113.1594),
    'M': AminoAcid('Methionine',     'M', 'Met', 131.04049, 131.1986),
    'N': AminoAcid('Asparagine',     'N', 'Asn', 114.04293, 114.1039),
    'O': AminoAcid('Pyrrolysine',    'O', 'Pyl', 255.15829, 255.3172),
    'P': AminoAcid('Proline',        'P', 'Pro',  97.05276,  97.1167),
    'Q': AminoAcid('Glutamine',      'Q', 'Gln', 128.05858, 128.1307),
    'R': AminoAcid('Arginine',       'R', 'Arg', 156.10111, 156.1875),
    'S': AminoAcid('Serine',         'S', 'Ser',  87.03203,  87.0782),
    'T': AminoAcid('Threonine',      'T', 'Thr', 101.04768, 101.1051),
    'U': AminoAcid('Selenocysteine', 'U', 'Sec', 150.95364, 150.0388),
    'V': AminoAcid('Valine',         'V', 'Val',  99.06841,  99.1326),
    'W': AminoAcid('Tryptophan',     'W', 'Trp', 186.07931, 186.2132),
    'Y': AminoAcid('Tyrosine',       'Y', 'Tyr', 163.06333, 163.1760)
}

skew_step={
    'A':0,
    'C':-1,
    'G': +1,
    'T': 0
}


 

# createSimpleDNASubst
#
# Populate a simple scoring table

# Inputs:   match    Reward for matching
#           subst    Penalty for a mismatch
#           bases    Replace with 'AUGC' for RNA
def createSimpleDNASubst(match=+1,subst=1,bases='ATGC'):
    weights={}
    for i in range(len(bases)):
        for j in range(len(bases)):
            weights[(bases[i],bases[j])] = +match if i==j else -subst          
    return weights

# get_re_protein
# Produce a regular expression to recognize a straing of amino acids

def get_re_protein(min_length=1):
    return re.compile('[A,C-IK-WY]{'+str(min_length)+',}')

if __name__=='__main__':
    import unittest
    
    class Test_Amino_acids(unittest.TestCase):
        def test_match_integer(self):
            for key in integer_masses:
                self.assertEqual(integer_masses[key],amino_acids[key].asInteger())
                
    unittest.main()