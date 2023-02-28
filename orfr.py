#!/usr/bin/env python
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
coding_dna = Seq("AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG", generic_dna)
print (coding_dna.translate())
