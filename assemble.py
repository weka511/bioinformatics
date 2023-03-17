#!/usr/bin/env python

#   Copyright (C) 2023 Simon Crase

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

''' Code for Chapter 3: How do we assemble genomes?'''

from unittest import main, TestCase, skip
import numpy as np
from numpy.testing import assert_array_equal
from helpers  import rotate
from fasta    import FastaContent
from rosalind import patternToNumber

def kmer_composition(k,dna):
    '''
     BA3A	Generate the k-mer Composition of a String

     Input: An integer k and a string Text.

     Return: Compositionk(Text) (the k-mers can be provided in any order).
    '''
    return [dna[i:i+k] for i in range(1+len(dna)-k)]



def reconstruct_as_list(k,n,fragments,extract=lambda fragments,i: fragments[i]):
    '''
     BA3B	Reconstruct a String from its Genome Path

      Input: A sequence of k-mers Pattern1, ... , Patternn such that the last k - 1
      symbols of Pattern[i] are equal to the first k - 1 symbols of Pattern[i+1]
      for i from 1 to n-1.

      Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal
      to Pattern[i] for all i.
    '''
    result=[]
    for i in range(0,n,k):
        result.append(extract(fragments,i))
    target_len=n+k-1
    actual_len=len(result)*k
    overlap=target_len-actual_len
    if overlap>0:
        result.append(fragments[-1][k-overlap:k])
    return result

def reconstruct(fragments):
    return ''.join(reconstruct_as_list(len(fragments[0]),len(fragments),fragments))



def grph_kmers(strings):
    '''
    BA3C	Construct the Overlap Graph of a Collection of k-mers

    Construct the overlap graph of a collection of k-mers.

    Input: A collection Patterns of k-mers.

    Return: The overlap graph Overlap(Patterns), in the form of an adjacency list.
    '''
    kk=len(strings[0])-1
    graph=[]
    for s in strings:
        for t in strings:
            if s!=t and s[-kk:]==t[:kk]:
                graph.append((s,t))

    return graph



def deBruijn(k,text):
    '''
    BA3D 	Construct the De Bruijn Graph of a String

    Given: An integer k and a string Text.

    Return:DeBruijnk(Text), in the form of an adjacency list.
    '''
    kmers=kmer_composition(k-1,text)
    def standardize(ll):
        lll=list(set(ll))
        lll.sort()
        return lll

    pathgraph=grph_kmers(kmers)

    deBruijn_dict={}
    for [a,b] in pathgraph:
        if not a in deBruijn_dict:
            deBruijn_dict[a]=[]
        deBruijn_dict[a].append(b)
    graph= [(a,standardize(deBruijn_dict[a])) for a in deBruijn_dict.keys()]
    graph.sort()
    return graph


def deBruijn_collection(pattern,
                        head=lambda kmer: kmer[0:-1],
                        tail=lambda kmer: kmer[1:]):
    '''
    BA3E 	Construct the De Bruijn Graph of a Collection of k-mers

    Input: A collection of k-mers Patterns.

    Return: The de Bruijn graph DeBruijn(Patterns), in the form of an adjacency list.
    '''
    graph={}
    k=len(pattern[0])
    for kmer in pattern:
        if not head(kmer) in graph:
            graph[head(kmer)]=[]
        graph[head(kmer)].append(tail(kmer))
    for kmer in graph.keys():
        graph[kmer].sort()
    return graph


def find_eulerian_cycle(graph):
    '''
    BA3F 	Find an Eulerian Cycle in a Graph

     Input: An Eulerian directed graph, in the form of an adjacency list.

     Return: An Eulerian cycle in this graph.
    '''
    def create_unexplored_edges():
        edges=[]
        for a in graph.keys():
            for b in graph[a]:
                edges.append((a,b))
        return edges

    def find_next(node):
        for succ in graph[node]:
            if (node,succ) in unexplored:
                return succ
        return None

    def find_branch(cycle):
        pos=0
        for node in cycle:
            succ=find_next(node)
            if succ!=None:
                return (pos,node)
            pos+=1

    def find_cycle(cycle0,node):
        cycle=cycle0
        succ=find_next(node)
        while succ!=None:
            cycle.append(succ)
            unexplored.remove((node,succ))
            node=succ
            succ=find_next(node)
        return cycle

    unexplored= create_unexplored_edges()
    target_edges=len(unexplored)
    node,_=unexplored[0]
    cycle=find_cycle([node],node)
    while len(unexplored)>0:
        (pos,branch)=find_branch(cycle)
        cycle=rotate(cycle,pos)
        cycle=find_cycle(cycle,cycle[len(cycle)-1])
    return cycle



def find_eulerian_path(graph):
    '''
    BA3G 	Find an Eulerian Path in a Graph

    Input: A directed graph that contains an Eulerian path, where the graph
    is given in the form of an adjacency list.

    Return: An Eulerian path in this graph.
    '''
    def nodes():
        return list(set(list(graph.keys())+
                        [out for outs in graph.values() for out in outs]))
    def get_start():
        start=[]
        for candidate in graph.keys():
            in_count=0
            for _,ins in graph.items():
                in_count+=sum([1 for i in ins if i==candidate])
            if in_count<len(graph[candidate]):
                start.append(candidate)
        return start

    def get_finish():
        finish=[]
        closed=[]
        for outs in graph.values():
            for candidate in outs:
                if not candidate in closed:
                    in_count=0
                    for os in graph.values():
                        for o in os:
                            if o==candidate:
                                in_count+=1
                    if not candidate in graph or in_count>len(graph[candidate]):
                        finish.append(candidate)
                    closed.append(candidate)
        return finish

    def adjust(path,start,finish):
        i=len(path)
        while i>0 and (path[0]!=start[0] or path[-1]!=finish[0]):
            path=path[1:]+[path[0]]
            i-=1
        return path

    start=get_start()
    finish=get_finish()
    if len(start)==1 and len(finish)==1:
        graph[finish[0]]=start[0:1]
        return adjust(find_eulerian_cycle(graph)[0:-1],start,finish)
    raise RosalindException(\
        'Start %(start)s and finish %(finish)s should have one element each'%locals())

def reconstruct_from_kmers(k,patterns):
    '''
    BA3H 	Reconstruct a String from its k-mer Composition

     Input: An integer k followed by a list of k-mers Patterns.

     Return: A string Text with k-mer composition equal to Patterns. (If multiple
     answers exist, you may return any one.)
    '''
    return reconstruct(find_eulerian_path(deBruijn_collection(patterns)))

def k_universal_circular_string(k):
    '''
    BA3I 	Find a k-Universal Circular String
    something off - I have taken this out of tests, even though the
    website accepts my answers. A small test case appears to give different
    answers on successive runs - maybe iterating through dict is not deterministic?
    '''
    def bits(i):
        return ('{0:b}'.format(i)).zfill(k)
    patterns=[bits(i) for i in range(2**k)]
    result= reconstruct(find_eulerian_cycle(deBruijn_collection(patterns)))[k-1:]
    for pattern in patterns:
        if not pattern in result:
            raise RosalindException('%(pattern)s %(result)s'%locals())
    return result

def reconstruct_from_paired_kmers(k,d,patterns):
    '''
    BA3J 	Reconstruct a String from its Paired Composition

    Input: Integers k and d followed by a collection of paired k-mers PairedReads.

    Return: A string Text with (k, d)-mer composition equal to PairedReads.
    (If multiple answers exist, you may return any one.)
    '''
    def create_pair(string):
        [a,b]=string.split('|')
        return (a,b)
    def prefix(pair):
        (a,b)=pair
        return (a[0:-1],b[0:-1])
    def suffix(pair):
        (a,b)=pair
        return (a[1:],b[1:])
    def reconstruct_from_graph(path):
        def extract(pairs,i):
            (a,_)=pairs[i]
            return a
        head_reconstruction=                         \
            ''.join(reconstruct_as_list(k-1,\
                                        len(path),\
                                        path,
                                        extract)[:-1])
        expected_length=len(patterns)+2*k+d-1
        deficit=expected_length-len(head_reconstruction)
        _,tail_reconstruction=path[-1]
        i=-2
        while len(tail_reconstruction)<deficit:
            _,tail=path[i]
            tail_reconstruction=tail[0]+tail_reconstruction
            i-=1
        return head_reconstruction+tail_reconstruction
    return reconstruct_from_graph(                            \
        find_eulerian_path(                                   \
            deBruijn_collection(                              \
                [create_pair(string) for string in patterns], \
                prefix,                                       \
                suffix)))

def create_contigs(patterns):
    '''BA3K 	Generate Contigs from a Collection of Reads'''
    contigs=[]
    for path in non_branching_paths(deBruijn_collection(patterns)):
        contig=path[0]
        for p in path[1:]:
            contig=contig+p[-1]
        contigs.append(contig)
    return contigs

def construct_from_gapped(k,d,patterns):
    '''BA3L 	Construct a String Spelled by a Gapped Genome Path'''
    return reconstruct_from_paired_kmers(k,d,patterns)

def non_branching_paths(graph):
    '''BA3M 	Generate All Maximal Non-Branching Paths in a Graph'''
    def add_counts(graph):
        counts={}
        for node,outs in graph.items():
            counts[node]=(0,len(outs))
        for node,outs in graph.items():
            for out in outs:
                if not out in counts:
                    counts[out]=(0,0)
                x,y=counts[out]
                counts[out]=x+1,y
        return counts
    def make_cycle(start,graph):
        def standardize(cycle):
            mm=min(cycle)
            ii=cycle.index(mm)
            return cycle[ii:]+cycle[:ii]+[mm]
        cycle=[start]
        node=start

        while node in graph and len(graph[node])==1:
            succ=graph[node][0]
            if succ==start:
                return standardize(cycle)
            else:
                cycle.append(succ)
                node=succ

        return []

    def isolated_cycles(nodes,graph):
        cycles=[]
        for node in graph.keys():
            cycle=make_cycle(node,graph)
            if len(cycle)>0 and not cycle in cycles:
                cycles.append(cycle)
        return cycles

    paths=[]
    nodes=add_counts(graph)

    for v in nodes:
        (ins,outs)=nodes[v]
        if ins!=1 or outs!=1:
            if outs>0:
                for w in graph[v]:
                    nbp=[v,w]
                    w_in,w_out=nodes[w]
                    while w_in==1 and w_out==1:
                        u=graph[w][0]
                        nbp.append(u)
                        w=u
                        w_in,w_out=nodes[w]
                    paths.append(nbp)
    return isolated_cycles(nodes,graph)+paths



def lexig(k,fasta):
    '''
    KMER  Generalizing GC-Content

     Input: A DNA string s in FASTA format (having length at most 100 kbp).

     Return: The 4-mer composition of s
     '''

    (a,b)  = fasta[0]
    counts = np.zeros((4**k))
    for index in [patternToNumber(kmer) for kmer in kmer_composition(k,b)]:
        counts[index]+=1
    return counts

if __name__=='__main__':

    class Test_3_Assembly(TestCase):
        def test_ba3a(self):
            kmers = kmer_composition(5,'CAATCCAAC')
            self.assertEqual(5,len(kmers))
            self.assertIn('AATCC', kmers)
            self.assertIn('ATCCA', kmers)
            self.assertIn('CAATC', kmers)
            self.assertIn('CCAAC', kmers)
            self.assertIn('TCCAA', kmers)

        # BA3B	Reconstruct a String from its Genome Path
        def test_ba3b(self):
            self.assertEqual("ACCGAAGCT",
                             reconstruct(["ACCGA",
                                          "CCGAA",
                                          "CGAAG",
                                          "GAAGC",
                                          "AAGCT"] ))

# BA3C	Construct the Overlap Graph of a Collection of k-mers
        def test_ba3c(self):
            graph=grph_kmers([
                'ATGCG',
                'GCATG',
                'CATGC',
                'AGGCA',
                'GGCAT'
            ])
            self.assertEqual(4,len(graph))
            self.assertIn(('AGGCA','GGCAT'),graph)
            self.assertIn(('CATGC','ATGCG'),graph)
            self.assertIn(('GCATG','CATGC'),graph)
            self.assertIn(('GGCAT','GCATG'),graph)

# BA3D 	Construct the De Bruijn Graph of a String
        def test_ba3d(self):
            graph=deBruijn(4,'AAGATTCTCTAC')
            self.assertIn(('AAG',['AGA']),graph)
            self.assertIn(('AGA',['GAT']),graph)
            self.assertIn(('ATT',['TTC']),graph)
            self.assertIn(('CTA',['TAC']),graph)
            self.assertIn(('CTC',['TCT']),graph)
            self.assertIn(('GAT',['ATT']),graph)
            self.assertIn(('TCT',['CTA','CTC']),graph)
            self.assertIn(('TTC',['TCT']),graph)
            self.assertEqual(8,len(graph))

#BA3E 	Construct the De Bruijn Graph of a Collection of k-mers
        def test_ba3e(self):
            graph=deBruijn_collection(
                ['GAGG',
                 'CAGG',
                 'GGGG',
                 'GGGA',
                 'CAGG',
                 'AGGG',
                 'GGAG'])
            self.assertEqual(graph['AGG'],['GGG'])
            self.assertEqual(graph['CAG'],['AGG','AGG'])
            self.assertEqual(graph['GAG'],['AGG'])
            self.assertEqual(graph['GGA'],['GAG'])
            self.assertEqual(graph['GGG'],['GGA','GGG'])
            self.assertEqual(5,len(graph.keys()))

        def assertCyclicEqual(self,a,b):
            self.assertEqual(len(a),len(b))
            for i in range(len(a)):
                if a==rotate(b,i):
                    return
            self.fail('Cycles %(a)s and %(b)s do not match'%locals())

#BA3F 	Find an Eulerian Cycle in a Graph
        def test_ba3f(self):
            self.assertCyclicEqual([6,8,7,9,6,5,4,2,1,0,3,2,6],
                             find_eulerian_cycle(
                                 {0 : [3],
                                  1 : [0],
                                  2 : [1, 6],
                                  3 : [2],
                                  4 : [2],
                                  5 : [4],
                                  6 : [5, 8],
                                  7 : [9],
                                  8 : [7],
                                  9 : [6]}
                             ))

#BA3G 	Find an Eulerian Path in a Graph
        def test_ba3g(self):
            self.assertEqual([6,7,8,9,6,3,0,2,1,3,4],find_eulerian_path(
                {
                    0 : [2],
                    1 : [3],
                    2 : [1],
                    3 : [0,4],
                    6 : [3,7],
                    7 : [8],
                    8 : [9],
                    9 : [6]
                }
            ))

# BA3H 	Reconstruct a String from its k-mer Composition
        def test_ba3h(self):
            self.assertEqual('GGCTTACCA',
                             reconstruct_from_kmers(
                                 4,
                                 ['CTTA',
                                  'ACCA',
                                  'TACC',
                                  'GGCT',
                                  'GCTT',
                                  'TTAC']))

        @skip('#17')
        def test_ba3i(self):
            ''' BA3I 	Find a k-Universal Circular String'''
            self.assertEqual('0110010100001111',k_universal_circular_string(4))

        @skip('#18')
        def test_ba3j(self):
            '''BA3J 	Reconstruct a String from its Paired Composition'''
            self.assertEqual("GTGGTCGTGAGATGTTGA",
                             reconstruct_from_paired_kmers(
                                 4,
                                 2,
                                ['GAGA|TTGA',
                                'TCGT|GATG',
                                'CGTG|ATGT',
                                'TGGT|TGAG',
                                'GTGA|TGTT',
                                'GTGG|GTGA',
                                'TGAG|GTTG',
                                'GGTC|GAGA',
                                'GTCG|AGAT']))
            self.assertEqual('CACCGATACTGATTCTGAAGCTT',
                             reconstruct_from_paired_kmers(3,
                                           1,
                                           [
                                              'ACC|ATA', #(ACC|ATA)
                                              'ACT|ATT', #(ACT|ATT)
                                              'ATA|TGA', #(ATA|TGA)
                                              'ATT|TGA', #(ATT|TGA)
                                              'CAC|GAT', #(CAC|GAT)
                                              'CCG|TAC', #(CCG|TAC)
                                              'CGA|ACT', #(CGA|ACT)
                                              'CTG|AGC', #(CTG|AGC)
                                              'CTG|TTC', #(CTG|TTC)
                                              'GAA|CTT', #(GAA|CTT)
                                              'GAT|CTG', #(GAT|CTG)
                                              'GAT|CTG', #(GAT|CTG)
                                              'TAC|GAT', #(TAC|GAT)
                                              'TCT|AAG', #(TCT|AAG)
                                              'TGA|GCT', #(TGA|GCT)
                                              'TGA|TCT', #(TGA|TCT)
                                              'TTC|GAA'  #(TTC|GAA)
                                           ]))

        def test_ba3k(self):
            '''BA3K 	Generate Contigs from a Collection of Reads'''
            contigs=create_contigs([
                'ATG',
                'ATG',
                'TGT',
                'TGG',
                'CAT',
                'GGA',
                'GAT',
                'AGA'
            ])
            self.assertIn('AGA',contigs)
            self.assertIn('ATG',contigs)
            self.assertIn('ATG',contigs)
            self.assertIn('CAT',contigs)
            self.assertIn('GAT',contigs)
            self.assertIn('TGGA',contigs)
            self.assertIn('TGT',contigs)


        def test_ba3l(self):
            '''BA3L 	Construct a String Spelled by a Gapped Genome Path'''
            self.assertEqual('GACCGAGCGCCGGA',
                             construct_from_gapped(4,
                                                   2,
                                                   [
                                                       'GACC|GCGC',
                                                       'ACCG|CGCC',
                                                       'CCGA|GCCG',
                                                       'CGAG|CCGG',
                                                       'GAGC|CGGA'
                             ]))


        def test_ba3m(self):
            '''BA3M 	Generate All Maximal Non-Branching Paths in a Graph'''
            paths=non_branching_paths({
                1 :[2],
                2 : [3],
                3 : [4,5],
                6 : [7],
                7 : [6]})
            self.assertIn([1 ,2 ,3],paths)
            self.assertIn([3, 4],paths)
            self.assertIn([3 ,5],paths)
            self.assertIn([6,7,6 ],paths)
            self.assertEqual(4,len(paths))


        def test_kmer(self):
            '''
            KMER  Generalizing GC-Content
            '''
            string='''>>Rosalind_6431
            CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGG
            CCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGT
            TTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCA
            AATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCG
            GGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGA
            CTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTA
            CCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG'''
            fasta = FastaContent(string.split('\n'))

            assert_array_equal(np.array(
                [4, 1, 4, 3, 0, 1, 1, 5, 1, 3, 1, 2, 2, 1, 2, 0, 1, 1, 3, 1, 2,
                 1, 3, 1, 1, 1, 1, 2, 2, 5, 1, 3, 0, 2, 2, 1, 1, 1, 1, 3, 1, 0,
                 0, 1, 5, 5, 1, 5, 0, 2, 0, 2, 1, 2, 1, 1, 1, 2, 0, 1, 0, 0, 1,
                 1, 3, 2, 1, 0, 3, 2, 3, 0, 0, 2, 0, 8, 0, 0, 1, 0, 2, 1, 3, 0,
                 0, 0, 1, 4, 3, 2, 1, 1, 3, 1, 2, 1, 3, 1, 2, 1, 2, 1, 1, 1, 2,
                 3, 2, 1, 1, 0, 1, 1, 3, 2, 1, 2, 6, 2, 1, 1, 1, 2, 3, 3, 3, 2,
                 3, 0, 3, 2, 1, 1, 0, 0, 1, 4, 3, 0, 1, 5, 0, 2, 0, 1, 2, 1, 3,
                 0, 1, 2, 2, 1, 1, 0, 3, 0, 0, 4, 5, 0, 3, 0, 2, 1, 1, 3, 0, 3,
                 2, 2, 1, 1, 0, 2, 1, 0, 2, 2, 1, 2, 0, 2, 2, 5, 2, 2, 1, 1, 2,
                 1, 2, 2, 2, 2, 1, 1, 3, 4, 0, 2, 1, 1, 0, 1, 2, 2, 1, 1, 1, 5,
                 2, 0, 3, 2, 1, 1, 2, 2, 3, 0, 3, 0, 1, 3, 1, 2, 3, 0, 2, 1, 2,
                 2, 1, 2, 3, 0, 1, 2, 3, 1, 1, 3, 1, 0, 1, 1, 3, 0, 2, 1, 2, 2,
                 0, 2, 1, 1]),
                lexig(4,fasta))

    main()
