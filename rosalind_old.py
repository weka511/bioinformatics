# Copyright (C) 2015-2019 Greenweaves Software Limited

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

import helpers as rh, reference_tables as rrt, fasta as f, math, sys, numpy,random, functools

class RosalindException(Exception):
    def __init__( self, message ):
        Exception.__init__(self, message)
        
## Rosalind functions



















### Which DNA Patterns play the role of molecular clocks

# BA2A 	Implement MotifEnumeration
#
# Input: Integers k and d, followed by a collection of strings Dna.
#
# Return: All (k, d)-motifs in Dna.

def enumerateMotifs(k,d,dna):
    approximate_matches={}
    kmers= rh.k_mers(k)
    
    def near_matches(kmer):
        if not kmer in approximate_matches:
            approximate_matches[kmer]=[kk for kk in kmers if hamm(kk,kmer)<=d]
        return approximate_matches[kmer]
    def good_enough(pattern):
        def match(string):
            for i in range(len(string)-k+1):
                kmer=string[i:i+k]
                if hamm(kmer,pattern)<=d:
                    return True
            return False
        for string in dna:
            if not match(string):
                return False
        return True
                    
    patterns=[]
    for string in dna:
        for i in range(len(string)-k+1):
            kmer=string[i:i+k]
            for pattern in near_matches(kmer):
                if good_enough(pattern):
                    patterns.append(pattern)
    return ' '.join(sorted(list(set(patterns))))

# BA2B 	Find a Median String 
#
# Input: An integer k and a collection of strings Dna.
#
# Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers
# Pattern. (If multiple answers exist, you may return any one.)

def medianString(k,dna):
    def findClosest(d): 
        distance=sys.float_info.max
        closest=None
        for k_mer in rh.k_mers(k):
            if distance>d(k_mer,dna):
                distance=d(k_mer,dna)
                closest=k_mer
        return closest
    return findClosest(distanceBetweenPatternAndStrings)


# BA2C 	Find a Profile-most Probable k-mer in a String	
#
# Input: A string Text, an integer k, and a 4 × k matrix Profile.
#
# Return: A Profile-most probable k-mer in Text. 

def mostProbable(text,n,profile):
    # probability of kmer given profile
    def prob(kmer):
        p=1
        for j in range(n):
            i=rrt.bases.find(kmer[j])
            p*=profile[i][j]
        return p
    
    def findMostProbable():
        probability=-1
        best=[]
        for i in range(len(text)-n+1):
            p=prob(text[i:i+n])
            if probability<p:
                probability=p
                best=text[i:i+n]
        return best
    
    return findMostProbable()

# BA2D 	Implement GreedyMotifSearch	
# BA2E 	Implement GreedyMotifSearch with Pseudocounts
#
# Input: Integers k and t, followed by a collection of strings Dna.
#        Optional parameter pseudo_counts specifies whether pseudo counts are to be used
#
# Return: A collection of strings BestMotifs resulting from running 
# GreedyMotifSearch(Dna, k, t). If at any step you find more than one
# Profile-most probable k-mer in a given string, use the one occurring first.
def greedyMotifSearch(k,t,dna,pseudo_counts=False):
    
    # Create an array containing the count of occurences of
    # each base at each position, summed over all motifs
    def count_occurrences_of_bases(
        motifs,\
        initialise_counts=numpy.ones if pseudo_counts else numpy.zeros
        ):
        matrix = initialise_counts((len(rrt.bases),k),dtype=int)
        for kmer in motifs:
            for j in range(k):
                i=rrt.bases.find(kmer[j])
                matrix[i,j]+=1
        return matrix
    
    def profile(motifs):
        return count_occurrences_of_bases(motifs)/float(len(motifs))
    
    def score(motifs):
        matrix=count_occurrences_of_bases(motifs)
        total=0
        for j in range(k):
            m=0
            for i in range(len(rrt.bases)):
                if m<matrix[i,j]:
                    m=matrix[i,j]
            total+=(len(rrt.bases)-m)
        return total
    
    bestMotifs=[genome[0:k] for genome in dna]
    for motif in [dna[0][i:i+k] for i in range(len(dna[0])-k+1)]:
        motifs=[motif]
        for i in range(1,t):
            motifs.append(mostProbable(dna[i],k,profile(motifs)))
        if score(motifs)<score(bestMotifs):
            bestMotifs=motifs
    return bestMotifs	  	  	 

	
	  	  	 
# BA2F 	Implement RandomizedMotifSearch	

def randomized_motif_search(k,t,dna,eps=1):
    def score(k,motifs):
        total=0
        for j in range(k):
            counts=numpy.zeros(len(rrt.bases),dtype=numpy.int32)
            for motif in motifs:
                i=rrt.bases.find(motif[j])
                counts[i]+=1
            max=-1
            ii=-1
            for i in range(len(rrt.bases)):
                if max<counts[i]:
                    ii=i
                    max=counts[ii]
            for i in range(len(rrt.bases)):
                if i!=ii:
                    total+=counts[i]
        return total    
    def counts(motifs):
        matrix=numpy.ones((len(rrt.bases),k),dtype=int)
        for i in range(len(rrt.bases)):
            for j in range(k):
                matrix[i,j]*=eps
        for kmer in motifs:
            for j in range(k):
                i=rrt.bases.find(kmer[j])
                matrix[i,j]+=1
        return matrix
    def Motifs(profile,dna):
        def get_k(Profile):
            return len(Profile[0])        
        def prob(kmer):
            p=1
            for j in range(k):
                i=rrt.bases.find(kmer[j])
                p*=profile[i][j]
            return p    
        k=get_k(profile)
        motifs=[]
        for s in dna:
            max_probability=-1
            most_probable_kmer=''
            for kmer in [s[i:i+k].upper() for i in range(len(s)-k+1)]:
                #if len(kmer)<k:
                    #print kmer,s,i
                if max_probability<prob(kmer):
                    max_probability=prob(kmer)
                    most_probable_kmer=kmer
            motifs.append(most_probable_kmer)
        return motifs    
    def Profile(motifs):
        matrix=counts(motifs)
        probabilities=numpy.zeros((len(rrt.bases),k),dtype=float)
        for i in range(len(rrt.bases)):
            for j in range(k):
                probabilities[i,j]=matrix[i,j]/float(len(motifs))
        return probabilities
  
    def random_kmer(string):
        i=random.randint(0,len(string)-k)
        return string[i:i+k]
 
    motifs=[]
    
    for i in range(t):
        motifs.append(random_kmer(dna[i]))
    bestMotifs=motifs
    while True:
        profile=Profile(motifs)
        motifs = Motifs(profile, dna)
        if score(k,motifs) < score(k,bestMotifs):
            bestMotifs = motifs
        else:
            return (score(k,bestMotifs),bestMotifs)

def randomized_motif_search_driver(k,t,dna,N=1000):
    best=sys.float_info.max
    mm=[]
    for i in range(N):
        (sc,motifs) =randomized_motif_search(k,t,dna)
        if sc<best:
            best=sc
            mm=motifs
        if i%100==0:
            print (i,best)
            for m in mm:
                print (m)
    return (best,mm)

	  	  	 
# BA2G 	Implement GibbsSampler	
#
# Input: Integers k, t, and N, followed by a collection of strings Dna.
#
# Return: The strings BestMotifs resulting from running 
# GibbsSampler(Dna, k, t, N) with 20 random starts.
# Remember to use pseudocounts!

def gibbs(k,t,n,dna,eps=1):
    def score(k,motifs):
        total=0
        for j in range(k):
            counts=numpy.zeros(len(rrt.bases),dtype=numpy.int32)
            for motif in motifs:
                i=rrt.bases.find(motif[j])
                counts[i]+=1
            max=-1
            ii=-1
            for i in range(len(rrt.bases)):
                if max<counts[i]:
                    ii=i
                    max=counts[ii]
            for i in range(len(rrt.bases)):
                if i!=ii:
                    total+=counts[i]
        return total    
    def random_kmer(string):
        i=random.randint(0,len(string)-k)
        return string[i:i+k]
    def dropOneMotif(motifs,i):
        return [motifs[j] for j in range(len(motifs)) if j!=i]
    def counts(motifs):
        matrix=numpy.ones((len(rrt.bases),k),dtype=int)
        for i in range(len(rrt.bases)):
            for j in range(k):
                matrix[i,j]*=eps
        for kmer in motifs:
            for j in range(k):
                i=rrt.bases.find(kmer[j])
                matrix[i,j]+=1
        return matrix    
    def Profile(motifs):
        matrix=counts(motifs)
        probabilities=numpy.zeros((len(rrt.bases),k),dtype=float)
        for i in range(len(rrt.bases)):
            for j in range(k):
                probabilities[i,j]=matrix[i,j]/float(len(motifs))
        return probabilities
    
    def probability(kmer,profile):
        p=1
        for j in range(len(kmer)):
            i=rrt.bases.find(kmer[j])
            p*=profile[i][j]
        return p
    
    def accumulate(probabilities):
        total=0
        cumulative=[]
        for p in probabilities:
            total+=p
            cumulative.append(total)
        return cumulative
    
    def generate(probabilities):
        accumulated=accumulate(probabilities)
        rr=accumulated[len(accumulated)-1]*random.random()
        i=0
        while accumulated[i]<=rr:
            i+=1
        return i
    
    motifs=[]
    
    for i in range(t):
        motifs.append(random_kmer(dna[i]))
    bestMotifs=motifs
    
    trace=[]
    best_score=sys.float_info.max
    for j in range(n):
        i=random.randint(0,t-1)
        profile=Profile(dropOneMotif(motifs,i))
        motif_index=generate([probability(dna[i][ll:ll+k],profile)\
                              for ll in range(len(dna[i])-k)])
        motifs[i]=dna[i][motif_index:motif_index+k]
        sc=score(k,motifs)
        if  sc< best_score:
            best_score=sc
            bestMotifs = motifs
        trace.append(best_score)
                   
    return (score(k,bestMotifs),bestMotifs,trace)

	  	  	 
# BA2H 	Implement DistanceBetweenPatternAndStrings 

def distanceBetweenPatternAndStrings (pattern,dna):
    # Extend Hamming distance to work with string of unequal length
    def hamming(pattern,genome):
        return min([hamm(pattern,genome[i:i+len(pattern)]) 
                    for i in range(len(genome)-len(pattern)+1)])        
    return sum([hamming(pattern,motif) for motif in dna])
    
 

# BA3B	Reconstruct a String from its Genome Path
#
# Input: A sequence of k-mers Pattern1, ... , Patternn such that the last k - 1
# symbols of Pattern[i] are equal to the first k - 1 symbols of Pattern[i+1]
# for i from 1 to n-1.
#
# Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal
# to Pattern[i] for all i.

def reconstruct_as_list(k,n,fragments,extract=lambda fragments,i: fragments[i]):
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

# BA3C	Construct the Overlap Graph of a Collection of k-mers
#
# Construct the overlap graph of a collection of k-mers.
#
# Input: A collection Patterns of k-mers.
#
# Return: The overlap graph Overlap(Patterns), in the form of an adjacency list.

def grph_kmers(strings):
    kk=len(strings[0])-1
    graph=[]
    for s in strings:
        for t in strings:
            if s!=t and s[-kk:]==t[:kk]:
                graph.append((s,t))
            
    return graph

# BA3D 	Construct the De Bruijn Graph of a String 
#
# Given: An integer k and a string Text.
#
# Return:DeBruijnk(Text), in the form of an adjacency list.

def deBruijn(k,text):
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

#BA3E 	Construct the De Bruijn Graph of a Collection of k-mers
#
# Input: A collection of k-mers Patterns.
#
# Return: The de Bruijn graph DeBruijn(Patterns), in the form of an adjacency list.

def deBruijn_collection(pattern,\
                        head=lambda kmer: kmer[0:-1],
                        tail=lambda kmer: kmer[1:]):
    graph={}
    k=len(pattern[0])
    for kmer in pattern:
        if not head(kmer) in graph:
            graph[head(kmer)]=[]
        graph[head(kmer)].append(tail(kmer))
    for kmer in graph.keys():
        graph[kmer].sort()
    return graph
 
# BA3F 	Find an Eulerian Cycle in a Graph 
#
# Input: An Eulerian directed graph, in the form of an adjacency list.
#
# Return: An Eulerian cycle in this graph.

def find_eulerian_cycle(graph):
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
        cycle=rh.rotate(cycle,pos)
        cycle=find_cycle(cycle,cycle[len(cycle)-1])
    return cycle

# BA3G 	Find an Eulerian Path in a Graph
#
# Input: A directed graph that contains an Eulerian path, where the graph
# is given in the form of an adjacency list.
#
# Return: An Eulerian path in this graph.

def find_eulerian_path(graph):
    def nodes():
        return list(set(list(graph.keys())+ \
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

# BA3H 	Reconstruct a String from its k-mer Composition 
#
# Input: An integer k followed by a list of k-mers Patterns.
#
# Return: A string Text with k-mer composition equal to Patterns. (If multiple
# answers exist, you may return any one.)

def reconstruct_from_kmers(k,patterns):
    return reconstruct(find_eulerian_path(deBruijn_collection(patterns)))

# BA3I 	Find a k-Universal Circular String 
# something off - I have taken this out of tests, even though the
# website accepts may answers. A small test case appears to give differnt
# answwrs on successive runs - maybe iterating through dict is not deterministic?

def k_universal_circular_string(k):
    def bits(i):
        return ('{0:b}'.format(i)).zfill(k)
    patterns=[bits(i) for i in range(2**k)]
    result= reconstruct(find_eulerian_cycle(deBruijn_collection(patterns)))[k-1:]
    for pattern in patterns:
        if not pattern in result:
            raise RosalindException('%(pattern)s %(result)s'%locals())   
    return result

# BA3J 	Reconstruct a String from its Paired Composition
#
# Input: Integers k and d followed by a collection of paired k-mers PairedReads.
#
# Return: A string Text with (k, d)-mer composition equal to PairedReads.
# (If multiple answers exist, you may return any one.)
def reconstruct_from_paired_kmers(k,d,patterns):
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

# BA3K 	Generate Contigs from a Collection of Reads 

def create_contigs(patterns):  
    contigs=[]
    for path in non_branching_paths(deBruijn_collection(patterns)):
        contig=path[0]
        for p in path[1:]:
            contig=contig+p[-1]
        contigs.append(contig)   
    return contigs

#BA3L 	Construct a String Spelled by a Gapped Genome Path
#   how does this differ from BA4J?

def construct_from_gapped(k,d,patterns):
    return reconstruct_from_paired_kmers(k,d,patterns)

#BA3M 	Generate All Maximal Non-Branching Paths in a Graph 

def non_branching_paths(graph):
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


### How do we sequence antibodies?



# BA4B	Find Substrings of a Genome Encoding a Given Amino Acid String
#
# There are three different ways to divide a DNA string into codons for
# translation, one starting at each of the first three starting positions of
# the string. These different ways of dividing a DNA string into codons are
# called reading frames. Since DNA is double-stranded, a genome has six reading
# frames (three on each strand).
#
# We say that a DNA string Pattern encodes an amino acid string Peptide if
# the RNA string transcribed from either Pattern or its reverse complement
# Pattern translates into Peptide.
#
# Input: A DNA string Text and an amino acid string Peptide.
#
# Return: All substrings of Text encoding Peptide (if any such substrings exist)

def findEncodings(text,peptide):
    def encodes(dna):
        try:
            return prot(dna_to_rna(dna))==peptide
        except KeyError:
            return False    
    candidates=[text[i:i+3*len(peptide)] for i in range(len(text)-3)]
    return [rna for rna in candidates if encodes(rna) or encodes(revc(rna))]
            
# BA4C	Generate the Theoretical Spectrum of a Cyclic Peptide
#
# The workhorse of peptide sequencing is the mass spectrometer, an expensive
# molecular scale that shatters molecules into pieces and then weighs the 
# resulting fragments. The mass spectrometer measures the mass of a molecule
# in daltons (Da); 1 Da is approximately equal to the mass of a single nuclear
# particle (i.e., a proton or neutron).
#
# We will approximate the mass of a molecule by simply adding the number of
# protons and neutrons found in the molecule’s constituent atoms, which yields
# the molecule’s integer mass. For example, the amino acid "Gly", which has 
# chemical formula C2H3ON, has an integer mass of 57,
# since 2·12 + 3·1 + 1·16 + 1·14 = 57. Yet 1 Da is not exactly equal to the mass
# of a proton/neutron, and we may need to account for different naturally
# occurring isotopes of each atom when weighing a molecule. As a result, 
# amino acids typically have non-integer masses (e.g., "Gly" has total mass
# equal to approximately 57.02 Da); for simplicity, however, we will work with
# the integer mass table given in Figure 1.
#
# The theoretical spectrum of a cyclic peptide Peptide, denoted
# Cyclospectrum(Peptide), is the collection of all of the masses of its
# subpeptides, in addition to the mass 0 and the mass of the entire peptide.
# We will assume that the theoretical spectrum can contain duplicate elements,
# as is the case for "NQEL", where "NQ" and "EL" have the same mass.
#
# Input: An amino acid string Peptide.
#
# Return: Cyclospectrum(Peptide).

def cycloSpectrum(peptide,mass=rrt.integer_masses):
    def get_pairs(index_range):
        n=len(index_range)
        return [(i,j) for i in index_range for j in range(i,i+n) if j!=i]
    augmented_peptide=peptide+peptide   # allows easy extraction of substrings
                                        # fromcyclic peptide
    spectrum=[rh.get_mass(augmented_peptide[a:b],mass)\
              for (a,b) in get_pairs(range(len(peptide)))]
    spectrum.append(rh.get_mass('',mass))
    spectrum.append(rh.get_mass(peptide,mass)) 
    spectrum.sort()
    return spectrum

# BA4D	Compute the Number of Peptides of Given Total Mass
#
# In Generate the Theoretical Spectrum of a Cyclic Peptide, we generated the
# theoretical spectrum of a known cyclic peptide. Although this task is
# relatively easy, our aim in mass spectrometry is to solve the reverse problem:
# we must reconstruct an unknown peptide from its experimental spectrum. 
# We will start by assuming that a biologist is lucky enough to generate an
# ideal experimental spectrum Spectrum, which is one coinciding with the
# peptide’s theoretical spectrum. Can we reconstruct a peptide whose
# theoretical spectrum is Spectrum?
#
# Denote the total mass of an amino acid string Peptide as Mass(Peptide).
# In mass spectrometry experiments, whereas the peptide that generated a 
# spectrum is unknown, the peptide’s mass is typically known and is denoted
# ParentMass(Spectrum). Of course, given an ideal experimental spectrum,
# Mass(Peptide) is given by the largest mass in the spectrum.
#
# A brute force approach to reconstructing a peptide from its theoretical
# spectrum would generate all possible peptides whose mass is equal to 
# ParentMass(Spectrum) and then check which of these peptides has theoretical
# spectra matching Spectrum. However, we should be concerned about the running
# time of such an approach: how many peptides are there having mass equal
# to ParentMass(Spectrum)?
#
# Input: An integer m.
#
# Return: The number of linear peptides having integer mass m.
#
# NB, treat peptide as a vector of masses, so amino acids with the same
# mass are the same

def count_peptides_linear(total_mass):
    cache=[]
    masses=list(set(rrt.integer_masses.values()))
    for target_mass in range(total_mass+1):
        total=0
        for amino_acid_mass in masses:
            residual_mass=target_mass-amino_acid_mass
            if residual_mass==0:
                total+=1
            elif residual_mass>0:
                total+=cache[residual_mass]
        cache.append(total)
        
    return total


def get_weight(peptide):
    return sum(rrt.amino_acids[amino_acid].mon_mass for amino_acid in peptide)

def expand(peptides,masses):
    return [peptide+[mass] for peptide in peptides for mass in masses]
 
def mass(peptide):
    return sum([weight for weight in peptide])

def parentMass(spectrum):
    return max(spectrum)

# BA4E 	Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum
#
# Input: A collection of (possibly repeated) integers Spectrum corresponding
#       to an ideal experimental spectrum.
#
# Return: An amino acid string Peptide such that Cyclospectrum(Peptide) = 
#        Spectrum (if such a string exists).

def find_cyclopeptide_sequence(spectrum):
        
    def consistent(peptide,spectrum):
        def count(element,spect):
            return len ([s for s in spect if s==element])
        peptide_spectrum=linearSpectrum(peptide)
        for element in peptide_spectrum:
            if count(element,peptide_spectrum)>count(element,spectrum):
                return False
        return True
    
    def linearSpectrum(peptide):
        def get_pairs():
            return [(i,j) for i in range(len(peptide)) for j in range(i+1,len(peptide))]
        
        result=[sum(peptide[a:b]) for (a,b) in get_pairs()]
        result.append(0)
        result.append(sum(peptide)) 
        result.sort()
        return result 
        
    def cycloSpectrum(peptide):
        def get_pairs(index_range):
            n=len(index_range)
            return [(i,j) for i in index_range for j in range(i,i+n) if j!=i]
        augmented_peptide=peptide+peptide
        result=[sum(augmented_peptide[a:b]) for (a,b) in get_pairs(range(len(peptide)))]
        result.append(0)
        result.append(sum(peptide)) 
        result.sort()
        return result  
      
    peptides=[[]]
    output=[]
    masses=list(set(rrt.integer_masses.values()))
    
    while len(peptides)>0:
        next_peptides=[]
        for peptide in expand(peptides,masses):
            if mass(peptide) == parentMass(spectrum):
                if cycloSpectrum(peptide) == spectrum:
                    output.append(peptide)
            else:
                if consistent(peptide,spectrum):
                    next_peptides.append(peptide)    
        peptides=next_peptides
    return output

# BA4F 	Compute the Score of a Cyclic Peptide Against a Spectrum 
def score(peptide,spectrum,spect_from_peptide=cycloSpectrum):
    return rh.countMatchesInSpectra(spect_from_peptide(peptide),spectrum)

# BA4G 	Implement LeaderboardCyclopeptideSequencing 
def leaderPeptide(n,                                            \
                  spectrum,                                     \
                  masses=list(set(rrt.integer_masses.values())),\
                  spect1=rh.linearSpectrum,                     \
                  spect2=rh.linearSpectrum):
   
    
    leaderPeptide=[]
    leaderBoard=[leaderPeptide]
    
    while len(leaderBoard)>0:
        newBoard=[]
        for peptide in expand(leaderBoard,masses):
            if mass(peptide)==parentMass(spectrum):
                if score(peptide,spectrum,spect2)>\
                   score(leaderPeptide,spectrum,spect2):
                    leaderPeptide=peptide
                newBoard.append(peptide)
            elif mass(peptide)>parentMass(spectrum):
                pass #peptide will be dropped from leader board
            else:
                newBoard.append(peptide)
        leaderBoard=trim(newBoard, spectrum, n,spect1)
    return leaderPeptide

# BA4H 	Generate the Convolution of a Spectrum 
def convolution (spectrum):
    def create_counts(diffs):
        counts={}
        for diff in diffs:
            if not diff in counts:
                counts[diff]=0
            counts[diff]+=1
        return counts        
    diffs=[abs(spectrum[i]-spectrum[j])  \
           for i in range(len(spectrum))  \
           for j in range(i+1,len(spectrum)) if spectrum[i]!=spectrum[j]]
    return sorted([(diff,count)                                      \
                   for diff,count in create_counts(diffs).items()],  \
                  key=lambda x: (-x[1], x[0]))

# BA4I 	Implement ConvolutionCyclopeptideSequencing
#
# Given: An integer M, an integer N, and a collection of
# (possibly repeated) integers Spectrum.
#
# Return: A cyclic peptide LeaderPeptide with amino acids taken only from the
# top M elements (and ties) of the convolution of Spectrum that fall between
# 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).
#
# NB: I had to sort spectrum to pass the testcase in the textbook.

def convolutionCyclopeptideSequencing(m,n,spectrum,low_mass=57,high_mass=200):
    def get_masses_from_spectrum():
        masses=[]
        last_count=0
        for mass,count in convolution(spectrum):
            if low_mass<=mass and mass<=high_mass:
                if len(masses)<m:
                    masses.append(mass)
                    last_count=count
                else:
                    if count==last_count:
                        masses.append(mass)
                    else:
                        return masses
        return masses
    
    return leaderPeptide(n,spectrum,get_masses_from_spectrum(),spect2=rh.cycloSpectrum1)

# BA4L 	Trim a Peptide Leaderboard	
#
# Input: A leaderboard of linear peptides Leaderboard, a linear spectrum
# Spectrum, and an integer N.
#
# Return: The top N peptides from Leaderboard scored against Spectrum.
# Remember to use LinearScore.

def trim(leaderBoard, spectrum,n,spectrum_generator=rh.linearSpectrum):
    if len(leaderBoard)<n:
        return leaderBoard
    peptides_with_scores=[\
        (score(peptide,spectrum,spectrum_generator),peptide)\
        for peptide in leaderBoard]
    peptides_with_scores.sort(reverse=True)   
    (cutoff,_)= peptides_with_scores[n-1]
    return [peptide                                     \
            for (score,peptide) in peptides_with_scores \
            if score>=cutoff]

# Adapter for 'trim', so it will work with peptides as strings

def trim_for_strings(leaderBoard, spectrum,n,\
                     spectrum_generator=rh.linearSpectrum,\
                     masses=rrt.integer_masses):
    numeric_peptides=[]
    trans={}
    for peptide in leaderBoard:
        num=[masses[amino_acid] for amino_acid in peptide]
        numeric_peptides.append(num)
        trans[tuple(num)]=peptide
    trimmed=trim(numeric_peptides, spectrum,n,spectrum_generator)
    return [trans [tuple(peptide)] for peptide in trimmed]

# BA4M 	Solve the Turnpike Problem 
def turnpike(differences):
    def diffs(seq):
        return sorted([x-y for x in seq for y in seq if x>y]) 
    def match(rs,ss):
        for r,s in zip(rs,ss):
            if r!=s:
                return False
        return True    
    def get_length():
        n_zeroes=0
        i=0
        j=len(differences)-1
        index_first_pos=-1
        while i<j:
            if differences[i]+differences[j]!=0:
                raise RosalindException('Mismatch %(i)d and %(j)d'%locals())
            if differences[i]==0:
                n_zeroes+=1
            else:
                index_first_pos=j
            i+=1
            j-=1
        n_zeroes = 2*n_zeroes if len(differences)%2==0 else 2*n_zeroes+1
        if n_zeroes*n_zeroes==len(differences):
            return (index_first_pos,n_zeroes)
        else:
            raise RosalindException('Mismatched lengths')
        

    index_first_pos,length_result =get_length()
    largest=differences[-1]
    indices=[]
    while len(indices)<length_result-2:
        indices.append(index_first_pos+1)
    
    while True:
        result=[0]
        for index in indices:
            result.append(differences[index])
        result.append(largest)
       
        if match(diffs(result),differences[index_first_pos:]):
            return result
        j=len(indices)-1
        while j>0:
            if indices[j]<len(differences)-1:
                indices[j]+=1
                break
            else:
                indices[j]=index_first_pos+1
                j-=1
            
        
# BA4J 	Generate the Theoretical Spectrum of a Linear Peptide 
def linearSpectrumFromString(peptide):
    return rh.linearSpectrum([rrt.integer_masses[a] for a in peptide])

# BA4K 	Compute the Score of a Linear Peptide 
def linearScore(peptide,spectrum):
    return rh.countMatchesInSpectra(linearSpectrumFromString(peptide),spectrum)

def convolution_expanded(spectrum):
    return [diff for (diff,count) in convolution (spectrum) for i in range(count) ]
	  	  	 

### Test cases ###

if __name__=='__main__':
 
    import unittest,fragile
    
    class TestRosalind(unittest.TestCase):
              
 
 

        
        def test_ba6e(self):
            pairs=fragile.find_shared_kmers(3,'AAACTCATC','TTTCAAATC')
            self.assertIn((0, 4),pairs)
            self.assertIn((0, 0),pairs)
            self.assertIn((4, 2),pairs)
            self.assertIn((6, 6),pairs)
    

        

 

   
  

# BA2A 	Implement MotifEnumeration

        def test_ba2a(self):
            self.assertEqual('ATA ATT GTT TTT',
                             enumerateMotifs(3,
                                   1,
                                   ['ATTTGGC',
                                    'TGCCTTA',
                                    'CGGTATC',
                                    'GAAAATT']))

# BA2B 	Find a Median String 

        def test_ba2b(self):
            self.assertEqual('GAC',
                             medianString(3,[
                                 'AAATTGACGCAT',
                                 'GACGACCACGTT',
                                 'CGTCAGCGCCTG',
                                 'GCTGAGCACCGG',
                                 'AGTACGGGACAG'                               
                             ]))

          
            
# BA2C 	Find a Profile-most Probable k-mer in a String
        def test_ba2c(self):
            self.assertEqual(
                'CCGAG',
                mostProbable(
                    'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT',\
                    5,
                    [[0.2, 0.2, 0.3, 0.2, 0.3],
                     [0.4, 0.3, 0.1, 0.5, 0.1],
                     [0.3, 0.3, 0.5, 0.2, 0.4],
                     [0.1, 0.2, 0.1, 0.1, 0.2]]))

# BA2D 	Implement GreedyMotifSearch
        def test_ba2d(self):
            motifs=greedyMotifSearch(3,
                                     5,
                                     [
                                         'GGCGTTCAGGCA',
                                         'AAGAATCAGTCA',
                                         'CAAGGAGTTCGC',
                                         'CACGTCAATCAC',
                                         'CAATAATATTCG'                
                                     ])
            self.assertEqual(['CAA','CAA','CAA','CAG','CAG'],sorted(motifs))

       
# BA2E 	Implement GreedyMotifSearch with Pseudocounts  
        def test_ba2e(self):
            motifs=greedyMotifSearch(3,
                                     5,
                                     [
                                         'GGCGTTCAGGCA',
                                         'AAGAATCAGTCA',
                                         'CAAGGAGTTCGC',
                                         'CACGTCAATCAC',
                                         'CAATAATATTCG'                
                                     ],
                                     pseudo_counts=True)
            self.assertEqual(['ATC','ATC','TTC','TTC','TTC'],sorted(motifs))

# BA2F 	Implement RandomizedMotifSearch	
        #def test_ba2f(self):
            #(c,x)=randomized_motif_search_driver(8, 5,[
                #'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                #'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                #'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                #'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                #'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'],100000)
            #print (c)
            #self.assertIn('TCTCGGGG',x)
            #self.assertIn('CCAAGGTG',x)
            #self.assertIn('TACAGGCG',x)
            #self.assertIn('TTCAGGTG',x)
            #self.assertIn('TCCACGTG',x)

# BA2G 	Implement GibbsSampler	
        def test_ba2g(self):
            x=gibbs(8, 5, 100,[
                'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'])

# BA2H 	Implement DistanceBetweenPatternAndStrings 
        def test_ba2h(self):
            self.assertEqual(
                5,
                distanceBetweenPatternAndStrings('AAA',
                                                 ['TTACCTTAAC',
                                                  'GATATCTGTC',
                                                  'ACGGCGTTCG',
                                                  'CCCTAAAGAG',
                                                  'CGTCAGAGGT']))
            
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
                if a==rh.rotate(b,i):
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

# BA3I 	Find a k-Universal Circular String 
        #def test_ba3i(self):
            #self.assertEqual('0110010100001111',k_universal_circular_string(4))

# BA3J 	Reconstruct a String from its Paired Composition 
        def test_ba3j(self):
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

# BA3K 	Generate Contigs from a Collection of Reads
        def test_ba3k(self):
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

#BA3L 	Construct a String Spelled by a Gapped Genome Path 
        def test_ba3l(self):
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
    
#BA3M 	Generate All Maximal Non-Branching Paths in a Graph 
        def test_ba3m(self):   
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

# BA4B	Find Substrings of a Genome Encoding a Given Amino Acid String
        def test_ba4b(self):
            encodings=findEncodings('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA','MA')
            self.assertEqual(3,len(encodings))
            self.assertIn('ATGGCC',encodings)
            self.assertIn('GGCCAT',encodings)
            self.assertIn('ATGGCC',encodings)

# BA4C	Generate the Theoretical Spectrum of a Cyclic Peptide
        def test_ba4c(self):
            self.assertEqual([0,113,114,128,129,227,242,242,257,355,356,370,371,484],
                             cycloSpectrum('LEQN'))
        
        def test_ba4d(self):
            self.assertEqual(14712706211,count_peptides_linear(1024))
        
        def test_chain(self):
            self.assertAlmostEqual(821.392,get_weight('SKADYEK'),places=3)
        

                
        def test_ba4e(self):
            seq=find_cyclopeptide_sequence([0,113,128,186,241,299,314,427])
            self.assertEqual(6,len(seq))
            self.assertIn([186,128,113],seq)
            self.assertIn([186,113,128],seq)
            self.assertIn([128,186,113],seq)
            self.assertIn([128,113,186],seq)
            self.assertIn([113,186,128],seq)
            self.assertIn([113,128,186],seq)
 
# BA4F 	Compute the Score of a Cyclic Peptide Against a Spectrum 
        def test_ba4f(self):
            self.assertEqual(
                11,
                score(
                    'NQEL',
                    [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]))

# BA4G 	Implement LeaderboardCyclopeptideSequencing 
        def test_ba4g(self):
            self.assertEqual([129, 71, 147, 113],
                             leaderPeptide(
                                 10,
                                 [0, 71, 113, 129, 147, 200, 218,\
                                  260, 313, 331, 347, 389, 460]))            

# BA4H 	Generate the Convolution of a Spectrum 
        def test_ba4h(self):
            self.assertEqual(
                [137, 137, 186, 186, 49, 323],
                convolution_expanded([0, 137, 186, 323]))
            
        #def test_random(self):
            #self.assertEqual([-5.737, -5.217, -5.263, -5.360, -5.958, -6.628, -7.009],
                             #random_genome('CGATACAA',[0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783]))


# BA4I 	Implement ConvolutionCyclopeptideSequencing 
        def test_ba4i(self):
            self.assertEqual([99,71,137,57,72,57],
                             convolutionCyclopeptideSequencing(
                                 20,
                                 60,
                                 [57, 57, 71, 99, 129, 137,\
                                  170, 186, 194, 208, 228,\
                                  265, 285, 299, 307, 323,\
                                  356, 364, 394, 422, 493])) 
            
# BA4J 	Generate the Theoretical Spectrum of a Linear Peptide 
        def test_ba4j(self):
            self.assertEqual(
                [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484],
                linearSpectrumFromString('NQEL'))     

# BA4L 	Trim a Peptide Leaderboard
        def test_ba4l(self):
            self.assertEqual(['LAST', 'ALST'],
                             trim_for_strings(['LAST', 'ALST', 'TLLT', 'TQAS'],
                                  [0,71,87,101,113,158,184,188,259,271,372],
                                  2))

# BA4M 	Solve the Turnpike Problem

        def test_ba4m(self):        
            self.assertEqual([0, 2, 4,7, 10],
                             turnpike([
                                 -10,-8, -7, -6, -5,
                                 -4, -3, -3, -2, -2,
                                 0, 0, 0, 0, 0, 
                                 2, 2, 3, 3, 4,
                                 5, 6, 7, 8, 10
                             ]))


# BA5G 	Compute the Edit Distance Between Two Strings 	 	
                         
# BA5H 	Find a Highest-Scoring Fitting Alignment of Two Strings 	 	
                         
# BA5I 	Find a Highest-Scoring Overlap Alignment of Two Strings 	 	
                         
# BA5J 	Align Two Strings Using Affine Gap Penalties 	 	
                         
# BA5K 	Find a Middle Edge in an Alignment Graph in Linear Space 	 	
                         
# BA5L 	Align Two Strings Using Linear Space 	 	
                         
# BA5M 	Find a Highest-Scoring Multiple Sequence Alignment 	 	
                         
# BA5N 	Find a Topological Ordering of a DAG 
        def test_ba5n(self):
            self.assertEqual([5, 4, 1, 2, 3],
                             topological_order({
                                 1 : [2],
                                 2 : [3],
                                 4 : [2],
                                 5 : [3]
                             }))
            
#LGIS 	Longest Increasing Subsequence
        def test_longestIncreasingSubsequence(self):
            (a,d)=longestIncreasingSubsequence(5,[5, 1, 4, 2, 3])
            self.assertEqual([1,2,3],a)
            self.assertEqual([5,4,3],d)
            (a,d)=longestIncreasingSubsequence(9,[8, 2, 1, 6, 5, 7, 4, 3, 9])
            self.assertEqual([1, 5, 7, 9],a)
            self.assertEqual([8, 6, 5, 4, 3],d)
        
        def test_orf(self):
            string='''>Rosalind_99
            AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'''
            fasta=f.FastaContent(string.split('\n'))            
            peptides=get_reading_frames(fasta)
            self.assertIn('MLLGSFRLIPKETLIQVAGSSPCNLS',peptides)
            self.assertIn('M',peptides)
            self.assertIn('MGMTPRLGLESLLE',peptides)
            self.assertIn('MTPRLGLESLLE',peptides)            
            self.assertEqual(4,len(peptides))
            #for peptide in get_reading_frames(f.FastaFile('./data/rosalind_orf(3).txt') ):
                #print (peptide)

        #def test_wfmd(self):
            #self.assertAlmostEqual(0.772, wfmd(4,6,2,1),3)
            
        #def test_q5(self):
            #for peptide in ['MTAI',
                            #'MLAT',
                            #'MAIT',
                            #'MTAL',
                            #'TALM',
                            #'TLAM']:
                #sp=cycloSpectrum(peptide)
                #if sp==[0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]:
                    #print (peptide, sp)
         
        #def test_q6(self):
            ## 0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333
            #spectrum=[0 ,71, 99, 101, 103, 128, 129, 199, 200, 204, 227, 230, 231,\
                #298, 303, 328, 330, 332, 333]
            #print (spectrum)
            #for peptide in ['TVQ',
                            #'CTV',
                            #'AVQ',
                            #'AQV',
                            #'TCE',
                            #'VAQ']:
                #masses=[rrt.integer_masses[a] for a in peptide]
                #masses.sort()
                #ls=rh.linearSpectrum(masses)
                #if rh.consistent(masses,spectrum):
                    #print (peptide)
                #else:
                    #print (peptide,masses,ls)
                    
    unittest.main(exit=False)
