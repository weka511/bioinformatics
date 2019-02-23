# $HeadURL: https://server/svn/sandbox/trunk/rosalind/reference_tables.py $
# $LastChangedDate: 2016-01-18 09:22:31 +1300 (Mon, 18 Jan 2016) $
# $LastChangedRevision: 762 $

# Copyright (C) 2015-2016 Greenweaves Software Pty Ltd

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
     'UAA': ';',   'CAA': 'Q',      'AAA': 'K',      'GAA': 'E',
     'UAG': ';',   'CAG': 'Q',      'AAG': 'K',      'GAG': 'E',
     'UGU': 'C',      'CGU': 'R',      'AGU': 'S',      'GGU': 'G',
     'UGC': 'C',      'CGC': 'R',      'AGC': 'S',      'GGC': 'G',
     'UGA': ';',   'CGA': 'R',      'AGA': 'R',      'GGA': 'G',
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

class AminoAcid:
    def __init__(self,name,short,abbrev,mon_mass,average_mass):
        self.name=name
        self.short=short
        self.abbrev=abbrev
        self.mon_mass=mon_mass
        self.average_mass=average_mass
    def __str__(self):
        return '%(name)s %(short)s %(abbrev)s %(int_mass)d %(mon_mass)f %(average_mass)f'%\
               {
                   'name':self.name,
                   'short' :self.short,
                   'abbrev' : self.abbrev,
                   'mon_mass' : self.mon_mass,
                   'average_mass' : self.average_mass,
                   'int_mass' : self.asInteger()
               }
    def asInteger(self):
        return int(self.mon_mass)
    
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

# BLOSUM2
#
# populate the BLOSUM62 scoring table as described at http://rosalind.info/problems/ba5e/

def createBLOSUM62():
    amino_acids=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',\
          'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    raw_weights=[
        [4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2],
        [0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2],
        [-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3],
        [-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2],
        [-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3],
        [0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3],
        [-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2],
        [-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1],
        [-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2],
        [-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1],
        [-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1],
        [-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2],
        [-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3],
        [-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1],
        [-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2],
        [1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2],
        [0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2],
        [0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1],
        [-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2],
        [-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7]
    ]
    weights={}
    for i in range(len(amino_acids)):
        for j in range(len(amino_acids)):
            weights[(amino_acids[i],amino_acids[j])]=raw_weights[i][j]            
    return weights
 
# PAM250
#
# populate the PAM250 scoring table
def createPAM250():
    amino_acids=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',\
          'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    raw_weights=[
        [ 2,  -2,  0,  0, -3,1,-1,-1,-1,-2,-1,0,1,0,-2,1,1,0,-6,-3],
        [-2,  12, -5, -5, -4,-3,-3,-2,-5,-6,-5,-4,-3,-5,-4,0,-2,-2,-8,0],
        [ 0,  -5,  4,  3, -6,1,1,-2,0,-4,-3,2,-1,2,-1,0,0,-2,-7,-4],
        [ 0,  -5,  3,  4, -5,0,1,-2,0,-3,-2,1,-1,2,-1,0,0,-2,-7,-4],
        [-3,  -4, -6, -5,  9,-5,-2,1,-5,2,0,-3,-5,-5,-4,-3,-3,-1,0,7],
        [ 1,  -3,  1,  0, -5,5,-2,-3,-2,-4,-3,0,0,-1,-3,1,0,-1,-7,-5],
        [-1,  -3,  1,  1, -2,-2,6,-2,0,-2,-2,2,0,3,2,-1,-1,-2,-3,0],
        [-1,  -2, -2, -2,  1,-3,-2,5,-2,2,2,-2,-2,-2,-2,-1,0,4,-5,-1],
        [-1,  -5,  0,  0, -5,-2,0,-2,5,-3,0,1,-1,1,3,0,0,-2,-3,-4],
        [-2,  -6, -4, -3,  2,-4,-2,2,-3,6,4,-3,-3,-2,-3,-3,-2,2,-2,-1],
        [-1,  -5, -3, -2,  0,-3,-2,2,0,4,6,-2,-2,-1,0,-2,-1,2,-4,-2],
        [ 0,  -4,  2,  1, -3,0,2,-2,1,-3,-2,2,0,1,0,1,0,-2,-4,-2],
        [ 1,  -3, -1, -1, -5,0,0,-2,-1,-3,-2,0,6,0,0,1,0,-1,-6,-5],
        [ 0,  -5,  2,  2, -5,-1,3,-2,1,-2,-1,1,0,4,1,-1,-1,-2,-5,-4],
        [-2,  -4, -1, -1, -4,-3,2,-2,3,-3,0,0,0,1,6,0,-1,-2,2,-4],
        [ 1,   0,  0,  0, -3,1,-1,-1,0,-3,-2,1,1,-1,0,2,1,-1,-2,-3],
        [ 1,  -2,  0,  0, -3,0,-1,0,0,-2,-1,0,0,-1,-1,1,3,0,-5,-3],
        [ 0,  -2,  -2, -2,-1,-1,-2,4,-2,2,2,-2,-1,-2,-2,-1,0,4,-6,-2],
        [-6,  -8,  -7, -7,0,-7,-3,-5,-3,-2,-4,-4,-6,-5,2,-2,-5,-6,17,0],
        [-3,   0,  -4, -4,7,-5,0,-1,-4,-1,-2,-2,-5,-4,-4,-3,-3,-2,0,10]       
    ]
    weights={}
    for i in range(len(amino_acids)):
        for j in range(len(amino_acids)):
            weights[(amino_acids[i],amino_acids[j])]=raw_weights[i][j]            
    return weights

# Populate a simple scoring table

def createSimpleDNASubst():
    bases=['A', 'T', 'G', 'C']
 
    weights={}
    for i in range(len(bases)):
        for j in range(len(bases)):
            weights[(bases[i],bases[j])] = +1 if i==j else -1          
    return weights

if __name__=='__main__':
    for a in amino_acids.values():
        print (a)
    for key in integer_masses:
        print (integer_masses[key],amino_acids[key].asInteger())
        if integer_masses[key]!=amino_acids[key].asInteger():
            print ("mismatch")
    print ('BLOSUM62')
    weights= createBLOSUM62()
    for pair in sorted(weights.keys()):
        print (pair,weights[pair])
    print ('PAM250')
    weights= createPAM250()
    for pair in sorted(weights.keys()):
        print (pair,weights[pair])    
    