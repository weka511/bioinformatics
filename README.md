# Bioinformatics

Code inspired by [Bioinformatics Algorithms: an Active Learning Approach](http://bioinformaticsalgorithms.com/) and 
from [Rosalind](http://rosalind.info). 
NB: functions generally use zero based indexing; [Rosalind](http://rosalind.info) uses 1-based.

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|1|-| Where in the Genome does DNA replication begin?|
|BA1A|rosalind.py|[Compute the Number of Times a Pattern Appears in a Text](http://rosalind.info/problems/ba1a/) |
|BA1B|rosalind.py|[Find the Most Frequent Words in a String](http://rosalind.info/problems/ba1b/) |
|BA1C|rosalind.py|[Find the Reverse Complement of a String](http://rosalind.info/problems/ba1c/) |
|BA1D|rosalind.py|[Find All Occurrences of a Pattern in a String](http://rosalind.info/problems/ba1d/) |
|BA1E|rosalind.py|[Find Patterns Forming Clumps in a String](http://rosalind.info/problems/ba1e/)|
|BA1F|rosalind.py|[Find a Position in a Genome Minimizing the Skew](http://rosalind.info/problems/ba1f/)|
|BA1G|rosalind.py|[Compute the Hamming Distance Between Two Strings](http://rosalind.info/problems/ba1g/)|
|BA1H|rosalind.py|[Find All Approximate Occurrences of a Pattern in a String](http://rosalind.info/problems/ba1h/) |
|BA1I|rosalind.py|[Find the Most Frequent Words with Mismatches in a String](http://rosalind.info/problems/ba1i/) |
|BA1J|rosalind.py|[Find Frequent Words with Mismatches and Reverse Complements](http://rosalind.info/problems/ba1j/) |
|BA1K|rosalind.py|[Generate the Frequency Array of a Strings](http://rosalind.info/problems/ba1f/) |
|BA1L|rosalind.py|[Implement PatternToNumber](http://rosalind.info/problems/ba1l/) |
|BA1M|rosalind.py|[Implement NumberToPattern](http://rosalind.info/problems/ba1m/) |
|BA1N|rosalind.py|[Generate the d-Neighborhood of a String](http://rosalind.info/problems/ba1n/) |
|-|count_kmers.py| Challenge Problems from [Bioinformatics Algorithms](http://bioinformaticsalgorithms.com/): determine number of times each kmer appears in DNA, allowing for reverse complements and d-neighbourhoods.|
|2|-|Which DNA Patterns Play the Role of Molecular Clocks?|
|BA2A|rosalind.py|[Implement MotifEnumeration](http://rosalind.info/problems/ba2a/) |
|BA2B |rosalind.py|[Find a Median String ](http://rosalind.info/problems/ba1j/)
|BA2C |rosalind.py|[Find a Profile-most Probable k-mer in a String	](http://rosalind.info/problems/ba2c/)
|BA2D |rosalind.py|[Implement GreedyMotifSearch](http://rosalind.info/problems/ba2d/)
|BA2E |rosalind.py|[Implement GreedyMotifSearch with Pseudocounts](http://rosalind.info/problems/ba2e/)
|BA2F |rosalind.py|[Implement RandomizedMotifSearch](http://rosalind.info/problems/ba2f/)	
|BA2G |rosalind.py|[Implement GibbsSampler](http://rosalind.info/problems/ba2g/)
|BA2H |rosalind.py|[Implement DistanceBetweenPatternAndStrings](http://rosalind.info/problems/ba2h/) 
|3|-|How Do We Assemble Genomes?|
|BA3A|rosalind.py|[Generate the k-mer Composition of a String](http://rosalind.info/problems/ba3a/) |
|BA3B|rosalind.py|[Reconstruct a String from its Genome Path](http://rosalind.info/problems/ba3a/)
|BA3C|rosalind.py|[Construct the Overlap Graph of a Collection of k-mers](http://rosalind.info/problems/ba3a/)
|BA3D |rosalind.py|[Construct the De Bruijn Graph of a String ](http://rosalind.info/problems/ba3a/)
|BA3E |rosalind.py|[Construct the De Bruijn Graph of a Collection of k-mers](http://rosalind.info/problems/ba3a/)
|BA3F |rosalind.py|[Find an Eulerian Cycle in a Graph](http://rosalind.info/problems/ba3a/)
|BA3G |rosalind.py|[Find an Eulerian Path in a Graph](http://rosalind.info/problems/ba3a/)
|BA3H |rosalind.py|[Reconstruct a String from its k-mer Composition](http://rosalind.info/problems/ba3a/)
|BA3I |rosalind.py|[Find a k-Universal Circular String](http://rosalind.info/problems/ba3a/)
|BA3J |rosalind.py|[Reconstruct a String from its Paired Composition](http://rosalind.info/problems/ba3a/)
|4|spectrum.py|How Do We Sequence Antibiotics?|
|BA4A|rosalind.py|[Translate an RNA String into an Amino Acid String](http://rosalind.info/problems/ba4a/) |
|BA4B|rosalind.py|[Find Substrings of a Genome Encoding a Given Amino Acid String](http://rosalind.info/problems/ba4b/) |
|BA4C|rosalind.py|[Generate the Theoretical Spectrum of a Cyclic Peptide](http://rosalind.info/problems/ba4c/) |
|BA4D|rosalind.py|[Compute the Number of Peptides of Given Total Mass](http://rosalind.info/problems/ba4d/) |
|BA4E|rosalind.py|[Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum](http://rosalind.info/problems/ba4e/) |
|BA4F|rosalind.py|[Compute the Score of a Cyclic Peptide Against a Spectrum](http://rosalind.info/problems/ba4f/) |
|BA4G|rosalind.py|[Implement LeaderboardCyclopeptideSequencing ](http://rosalind.info/problems/ba4g/) |
|BA4H|rosalind.py|[Generate the Convolution of a Spectrum ](http://rosalind.info/problems/ba4h/) |
|BA4I|rosalind.py|[Implement ConvolutionCyclopeptideSequencing](http://rosalind.info/problems/ba4i/) |
|BA4J|rosalind.py|[Generate the Theoretical Spectrum of a Linear Peptide](http://rosalind.info/problems/ba4j/) |
|BA4K|rosalind.py|[Compute the Score of a Linear Peptide ](http://rosalind.info/problems/ba4k/) |
|BA4L|rosalind.py|[Trim a Peptide Leaderboard](http://rosalind.info/problems/ba4l/) |
|BA4M|ba4m.py spectrum.py| [Solve the Turnpike Problem](http://rosalind.info/problems/ba4m/) |
|-|tyrocidine.py| Challenge Problems from [Bioinformatics Algorithms](http://bioinformaticsalgorithms.com/): sequence Tyrocidine from supplied spectrum (WIP)|
|5|align.py|How do we compare Biological Sequences?|
|BA5A|align.py|[Find the Minimum Number of Coins Needed to Make Change](http://rosalind.info/problems/ba5a/)|
|BA5B|align.py|[Find the Length of a Longest Path in a Manhattan-like Grid](http://rosalind.info/problems/ba5b/)|
|BA5C|align.py|[Find a Longest Common Subsequence of Two Strings](http://rosalind.info/problems/ba5c/)|
|BA5D|align.py|[Find the Longest Path in a DAG](http://rosalind.info/problems/ba5d/)|
|BA5E|align.py|[Find a Highest-Scoring Alignment of Two Strings](http://rosalind.info/problems/ba5e/)|
|BA5F|align.py|[Find a Highest-Scoring Local Alignment of Two Strings](http://rosalind.info/problems/ba5f/)|
|BA5G|align.py| [Compute the Edit Distance Between Two Strings](http://rosalind.info/problems/ba5g/) |
|BA5H|BA5H.py align.py| [Find a Highest-Scoring Fitting Alignment of Two Strings](http://rosalind.info/problems/ba5h/) |   
|BA5I|BA5I.py align.py|[Find a Highest-Scoring Overlap Alignment of Two Strings](http://rosalind.info/problems/ba5i/)|    
|BA5J|BA5J.py align.py|[Align Two Strings Using Affine Gap Penalties](http://rosalind.info/problems/ba5j/)|    
|BA5K|BA5K.py align.py|[Find a Middle Edge in an Alignment Graph in Linear Space](http://rosalind.info/problems/ba5k/) |
|BA5L|BA5L.py align.py|[Align Two Strings Using Linear Space](http://rosalind.info/problems/ba5l/) (WIP)|    
|BA5M|BA5M.py align.py|[Find a Highest-Scoring Multiple Sequence Alignment](http://rosalind.info/problems/ba5m/) | 
|BA5N |align.py|[Find a Topological Ordering of a DAG](http://rosalind.info/problems/ba5n/) | 
|6|fragile.py|Are there fragile regions in the human genome?|
|BA6A|BA6A.py fragile.py|[Implement GreedySorting to Sort a Permutation by Reversals ](http://rosalind.info/problems/ba6a/)|
|BA6B|BA6B.py fragile.py|[Compute the Number of Breakpoints in a Permutation](http://rosalind.info/problems/ba6b/)|
|BA6C|BA6C.py fragile.py|[Compute the 2-Break Distance Between a Pair of Genomes](http://rosalind.info/problems/ba6c/)|
|BA6D|BA6D.py|[Find a Shortest Transformation of One Genome into Another by 2-Breaks](http://rosalind.info/problems/ba6d/) (wip)|
|BA6E|fragile.py| [Find All Shared k-mers of a Pair of Strings ](http://rosalind.info/problems/ba6e/)|
|BA6F|BA6F.py fragile.py|[Implement Chromosome to Cycle](http://rosalind.info/problems/ba6f/)|
|BA6G|BA6G.py fragile.py|[Implement Cycle to Chromosome](http://rosalind.info/problems/ba6g/)|
|BA6H|BA6H.py fragile.py|[Implement Coloured Edges](http://rosalind.info/problems/ba6h/)|
|BA6I|BA6I.py fragile.py|[Implement GraphToGenome ](http://rosalind.info/problems/ba6i/)|
|BA6J|BA6J.py fragile.py|[Implement 2-BreakOnGenomeGraph ](http://rosalind.info/problems/ba6j/)|
|BA6K|BA6K.py fragile.py|[Implement 2-BreakOnGenome ](http://rosalind.info/problems/ba6k/)|
|7|-|Which animal gave us SARS?|
|BA7A|BA7A.py|[Compute Distances Between Leaves](http://rosalind.info/problems/ba7a/)|
|BA7B|BA7B.py|[Limb Length Problem](http://rosalind.info/problems/ba7b/)|
|BA7C|BA7C.py|[Implement Additive Phylogeny](http://rosalind.info/problems/ba7c/)|
|BA7D|BA7D.py|[Implement UPGMA](http://rosalind.info/problems/ba7d/)|
|BA7E|BA7E.py|[Implement the Neighbor Joining Algorithm](http://rosalind.info/problems/ba7e/)|
|BA7F|BA7F.py|[Implement SmallParsimony](http://rosalind.info/problems/ba7f/)|
|BA7G|BA7G.py|[Adapt SmallParsimony to Unrooted Trees](http://rosalind.info/problems/ba7g/)|
|8|-|How did Yeast become a Winemaker?|
|BA8A|BA8A.py | [Implement FarthestFirstTraversal](http://rosalind.info/problems/ba8a/) |
|BA8B|BA8B.py | [Compute the Squared Error Distortion](http://rosalind.info/problems/ba8b/) |
|BA8C|BA8C.py | [Implement the Lloyd Algorithm for k-Means Clustering](http://rosalind.info/problems/ba8c/) |
|BA8D|BA8D.py | [Implement the Soft k-Means Clustering Algorithm](http://rosalind.info/problems/ba8d/) |
|BA8E|BA8E.py | [Implement Hierarchical Clustering](http://rosalind.info/problems/ba8e/) |
|9|snp.py|How do we locate disease causing mutations?|
|BA9A|BA9A.py snp.py| [Construct a Trie from a Collection of Patterns](http://rosalind.info/problems/ba9a/) |
|BA9B|BA9B.py snp.py| [Implement TrieMatching](http://rosalind.info/problems/ba9b/) |
|BA9C|BA9C.py snp.py| [Construct the Suffix Tree of a String](http://rosalind.info/problems/ba9c/)|
|BA9D|BA9D.py| [Find the Longest Repeat in a String](http://rosalind.info/problems/ba9d/)|
|BA9E|BA9E.py| [Find the Longest Substring Shared by Two Strings](http://rosalind.info/problems/ba9e/) (WIP)|
|BA9F|BA9F.py| [Find the Shortest Non-Shared Substring of Two Strings](http://rosalind.info/problems/ba9f/) (WIP)|
|BA9I|BA9I.py snp.py| [Construct the Burrows-Wheeler Transform of a String](http://rosalind.info/problems/ba9i/) |
|BA9J|BA9J.py snp.py| [Reconstruct a String from its Burrows-Wheeler Transform](http://rosalind.info/problems/ba9j/) |
|BA9K|BA9K.py snp.py| [Generate the Last-to-First Mapping of a String](http://rosalind.info/problems/ba9k/) (WIP)|
|10|hmm.py|Why have biologists still not developed an HIV vaccine?|
|BA10A|BA10A.py hmm.py| [Compute the Probability of a Hidden Path](http://rosalind.info/problems/ba10a/) |
|BA10B|BA10B.py hmm.py| [Compute the Probability of an Outcome Given a Hidden Path](http://rosalind.info/problems/ba10b/) |
|BA10C|BA10C.py hmm.py| [Implement the Viterbi Algorithm](http://rosalind.info/problems/ba10c/) |
|11|spectrum.py|Was T.rex just a big chicken?|
|BA11A|BA11A.py spectrum.py| [Spectrum Graph Construction](http://rosalind.info/problems/ba11a/) |
|BA11B|BA11B.py spectrum.py| [Implement DecodingIdealSpectrum](http://rosalind.info/problems/ba11b/) |
|BA11C|BA11C.py spectrum.py| [Convert a Peptide into a Peptide Vector ](http://rosalind.info/problems/ba11c/) |
|BA11D|BA11D.py spectrum.py| [Convert a Peptide Vector into a Peptide ](http://rosalind.info/problems/ba11d/) |
|fibo|fibo.py|[Fibonacci numbers](http://rosalind.info/problems/fibo/) |
|-|align.py|[Alignment](http://rosalind.info/problems/topics/alignment/)|
|edit|edit.py|[Edit Distance](http://rosalind.info/problems/edit/) |
|edta|edta.py|[Edit Distance Alignment](http://rosalind.info/problems/edta/)|
|gaff|gaff.py|[Global Alignment with Scoring Matrix and Affine Gap Penalty](http://rosalind.info/problems/gaff/)|
|gcon|gcon.py|[Global Alignment with Constant Gap Penalty](http://rosalind.info/problems/gcon/)|
|glob|GLOB.py|[Global Alignment with Scoring Matrix](http://rosalind.info/problems/glob/) -- test harness only|
|loca|loca.py|[Local Alignment with Scoring Matrix ](http://rosalind.info/problems/loca/) |
|hamm|rosalind.py|[Counting Point Mutations](http://rosalind.info/problems/hamm/) |
|oap|oap.py|[Overlap Alignment](http://rosalind.info/problems/oap/) (wip)|
|pdst|rosalind.py|[Creating a Distance Matrix ](http://rosalind.info/problems/pdst/) |
|sims|sims.py| [Finding Mutated Motifs](http://rosalind.info/problems/sims/)|
|tran|rosalind.py|[Transitions and Transversions](http://rosalind.info/problems/tran/) |
|-|-|[Combinatorics](http://rosalind.info/problems/topics/combinatorics/)|
|aspc|rosalind.py|[Introduction to Alternative Splicing](http://rosalind.info/problems/aspc/) |
|cat|rosalind.py|[Catalan Numbers and RNA Secondary Structures ](http://rosalind.info/problems/cat/) (WIP)|
|fibd|rosalind.py|[Mortal Fibonacci Rabbits](http://rosalind.info/problems/fibd/) |
|fib|rosalind.py|[Rabbits and Recurrence Relations](http://rosalind.info/problems/fib/) |
|mrna|rosalind.py|[Inferring mRNA from Protein](http://rosalind.info/problems/mrna/) |
|orf|orf.py|[Open Reading Frames](http://rosalind.info/problems/orf/)|
|pmch|pmch.py|[Perfect Matchings and RNA Secondary Structures](http://rosalind.info/problems/pmch/) |
|pper|rosalind.py|[Partial Permutations](http://rosalind.info/problems/pper/) |
|-|-|[Genome Assembly](http://rosalind.info/problems/topics/genome-assembly/)|
|corr|CORR.py|[Error Correction in Reads](http://rosalind.info/problems/corr/) |
|asmq|ASMQ.py|[Assessing Assembly Quality with N50 and N75](http://rosalind.info/problems/asmq/) |
|gasm|GASM.py|[Genome Assembly Using Reads](http://rosalind.info/problems/gasm/) |
|grep|GREP.py|[Genome Assembly with Perfect Coverage and Repeats](http://rosalind.info/problems/grep/)|
|long|LONG.py|[Genome Assembly as Shortest Superstring](http://rosalind.info/problems/long/) |
|pcov|PCOV.py|[Genome Assembly with Perfect Coverage](http://rosalind.info/problems/pcov/) |
|-|-|[Genome Rearrangements](http://rosalind.info/problems/topics/genome-rearrangements/)|
|lgis|lgis.py|[Longest Increasing Subsequence](http://rosalind.info/problems/lgis/) |
|rear|rear.py fragile.py|[Reversal distance](http://rosalind.info/problems/rear/) |
|sort|sort.py fragile.py|[Sorting by Reversals](http://rosalind.info/problems/sort/) |
|-|graphs.py|[Graphs](http://rosalind.info/problems/topics/graphs/) & [Graph Algorithms](http://rosalind.info/problems/topics/graph-algorithms/)|
|2sat|sat.py graphs.py|[2-Satisfiability](http://rosalind.info/problems/2sat/)|
|bf|bf.py graphs.py|[Bellman-Ford Algorithm](http://rosalind.info/problems/bf/)|
|bfs|bfs.py|[Beadth First search](http://rosalind.info/problems/bfs/) |
|bins|BINS.py|[Binary Search](http://rosalind.info/problems/bfs/) |
|bip|bip.py graphs.py|[Testing Bipartiteness](http://rosalind.info/problems/bip/) |
|cc|cc.py graphs.py|[Connected Components](http://rosalind.info/problems/cc/) |
|cte|cte.py graphs.py|[Shortest Cycle Through a Given Edge](http://rosalind.info/problems/cte/)|
|dag|dag.py graphs.py|[Testing Acyclicity](http://rosalind.info/problems/dag/) |
|ddeg|DDEG.py graphs.py|[Double-Degree Array](http://rosalind.info/problems/ddeg/) |
|deg|DEG.py graphs.py|[Degree Array](http://rosalind.info/problems/deg/) |
|dfs| graphs.py|[Depth First Search](http://rosalind.info/problems/dfs/) ???|
|dij|dij.py graphs.py|[Dijkstra's Algorithm](http://rosalind.info/problems/dij/) |
|grph|graphs.py|[Overlap Graphs](http://rosalind.info/problems/grph/) |
|gs|gs.py graphs.py|[General Sink](http://rosalind.info/problems/gs/)|
|hdag|hdag.py graphs.py|[Hamiltonian Path in DAG](http://rosalind.info/problems/hdag/) |
|nwc|nwc.py graphs.py|[Negative Weight Cycle](http://rosalind.info/problems/nwc/)|
|scc|scc.py graphs.py|[Strongly Connected Components](http://rosalind.info/problems/scc/)|
|sc|sc.py graphs.py|[Semi-Connected Graph](http://rosalind.info/problems/sc/)|
|sdag|sdag.py graphs.py|[Shortest Paths in DAG](http://rosalind.info/problems/sdag/)|
|sq|sq.py graphics.py|[Square in a Graph](http://rosalind.info/problems/sq/)|
|trie|rosalind.py|[Introduction to Pattern Matching](http://rosalind.info/problems/trie/) |
|ts|ts.py align.py|[Topological sort](http://rosalind.info/problems/ts/)|
|-|-|[Heredity](http://rosalind.info/problems/topics/heredity/)|
|iev|rosalind.py| [Calculating Expected Offspring](http://rosalind.info/problems/iev/)|
|indq|rosalind.py|[Independent Segregation of Chromosomes](http://rosalind.info/problems/indq/) |
|iprb|rosalind.py|[Mendel's First Law](http://rosalind.info/problems/iprb/) |
|lia|rosalind.py|[Independent Alleles](http://rosalind.info/problems/lia/) |
|sexl|rosalind.py|[Sex-Linked Inheritance ](http://rosalind.info/problems/sexl/) |
|-|phylogeny.py|[Phylogeny](http://rosalind.info/problems/topics/phylogeny/)|
|chbp|chbp.py phylogeny.py|[Character-Based Phylogeny](http://rosalind.info/problems/chbp/) WIP|
|cstr|cstr.py phylogeny.py|[Creating a Character Table from Genetic Strings](http://rosalind.info/problems/cstr/)|
|ctbl|ctbl.py phylogeny.py|[Creating a Character Table](http://rosalind.info/problems/ctbl/) |
|mend|mend.py|[Inferring a Genotype from a Pedigree](http://rosalind.info/problems/mend/) |
|nwck|nwck.py|[Distances in Trees](http://rosalind.info/problems/nwck/) |
|tree|tree.py|[Completing a Tree ](http://rosalind.info/problems/tree/) |
|-|-|[Population Dynamics](http://rosalind.info/problems/topics/population-dynamics/)|
|afrq|rosalind.py|[Counting Disease Carriers](http://rosalind.info/problems/afrq/) |
|foun|rosalind.py|[The Founder Effect and Genetic Drift](http://rosalind.info/problems/foun/) |
|wfmd|rosalind.py|[The Wright-Fisher Model of Genetic Drift](http://rosalind.info/problems/wfmd/) |
|-|-|[Probability](http://rosalind.info/problems/topics/probability/)|
|ebin|rosalind.py|[Wright-Fisher's Expected Behavior](http://rosalind.info/problems/ebin/) |
|eval|rosalind.py|[Expected Number of Restriction Sites](http://rosalind.info/problems/eval/) |
|rstr|rosalind.py|[Matching Random Motifs](http://rosalind.info/problems/rstr/) |
|-|-|Sequence Analysis|
|frmt|FRMT.py|[Data Formats](http://rosalind.info/problems/frmt/)|
|gbk|GBK.py|[GenBank Introduction](http://rosalind.info/problems/gbk/) |
|orfr|orfr.py|[Finding Genes with Open Reading Frames](http://rosalind.info/problems/orfr/) |
|ptra|PTRA.py|[Protein Translation](http://rosalind.info/problems/ptra/) |
|-|-|[Sorting](http://rosalind.info/problems/topics/sorting/)|
|2sum|2sum.py|[2SUM](http://rosalind.info/problems/2sum/) |
|3sum|3sum.py|[3SUM](http://rosalind.info/problems/3sum/) |
|ins|INS.py|[Insertion Sort](http://rosalind.info/problems/ins/) |
|inv|INV.py|[Counting Inversions](http://rosalind.info/problems/inv/) |
|med|MED.py|[Median](http://rosalind.info/problems/med/)|
|mer|MER.py|[Merge Two Sorted Arrays](http://rosalind.info/problems/mer/) |
|mprt|MPRT.py|[Finding a Protein Motif](http://rosalind.info/problems/mprt/)|
|ms|MS.py|[Mergesort](http://rosalind.info/problems/ms/)|
|par|PAR.py|[2-Way Partition](http://rosalind.info/problems/par/) |
|par3|PAR3.py|[3-Way Partition](http://rosalind.info/problems/par3/) |
|ps|ps.py|[Partial sort](http://rosalind.info/problems/ps/)|
|qs|QS.py|[Quicksort](http://rosalind.info/problems/qs/)|
|-|-|Strings & [String Algorithms](http://rosalind.info/problems/topics/string-algorithms/)|
|cons|rosalind.py|[Consensus and Profile](http://rosalind.info/problems/cons/) |
|dna|rosalind.py|[Counting DNA Nucleotides](http://rosalind.info/problems/dna/) |
|gc|rosalind.py|[Computing GC Content](http://rosalind.info/problems/gc/) |
|kmer|rosalind.py|[Generalizing GC-Content](http://rosalind.info/problems/kmer/) |
|kmp|kmp.py|[Speeding Up Motif Finding](http://rosalind.info/problems/kmp/) |
|lcsm|rosalind.py|[Finding a Shared Motif](http://rosalind.info/problems/lcsm/) |
|lcsq|lcsq.py|[Finding a Shared Spliced Motif](http://rosalind.info/problems/lcsq/) |
|ling|ling.py|[Linguistic Complexity of a Genome](http://rosalind.info/problems/ling/) (WIP)|
|revc|rosalind.py| [Complementing a Strand of DNA](http://rosalind.info/problems/revc/)|
|prot|rosalind.py|[Translating RNA into Protein](http://rosalind.info/problems/prot/) |
|revp|rosalind.py| [Locating Restriction Sites](http://rosalind.info/problems/revp/)|
|rna|rosalind.py|[Transcribing DNA into RNA](http://rosalind.info/problems/rna/) |
|scsp|scsp.py|[Interleaving Two Motifs](http://rosalind.info/problems/scsp/)|
|splc|rosalind.py|[RNA Splicing](http://rosalind.info/problems/splc/) |
|sseq|rosalind.py|[Finding a Spliced Motif](http://rosalind.info/problems/sseq/) |
|subs|rosalind.py|[Finding a Motif in DNA](http://rosalind.info/problems/subs/) |
|-|-|[Computational Mass Spectrometry](http://rosalind.info/problems/topics/computational-mass-spectrometry/)|
|conv|conv.py+spectrum.py|[Comparing Spectra with the Spectral Convolution](http://rosalind.info/problems/conv/)|
|full|full.py+spectrum.py|[Inferring Peptide from Full Spectrum](http://rosalind.info/problems/full/)|
|prtm|rosalind.py|[Calculating Protein Mass](http://rosalind.info/problems/prtm/) - get_weight(...)|
|prsm|prsm.py+spectrum.py|[Matching a Spectrum to a Protein](http://rosalind.info/problems/prsm/)|
|spec|spec.py+spectrum.py|[Inferring Protein from Spectrum](http://rosalind.info/problems/spec/)|
|sgra|sgra.py+spectrum.py|[Using the Spectrum Graph to Infer Peptides](http://rosalind.info/problems/sgra/)|
|-|-|[Set Theory](http://rosalind.info/problems/topics/set-theory/)|
|pdpl|pdpl|[Creating a Restriction Map](http://rosalind.info/problems/pfpl/)|
|seto|-|[Set operations](http://rosalind.info/problems/seto/) Solved using Python shell|
|sset|-|[Counting Subsets](http://rosalind.info/problems/sset/) Solved with a single line in Python shell: 2**n%1000000|
|-|-|Support| 
|-|helpers.py|Utilities for formatting output, parsing input, etc |
|-|LICENSE|License Agreement|
|-|newick.py|Parser for files in [Newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html) |
|-|README.md|This file|
|-|rosalind.py|Shared code|
|-|rosalind.wpr|WingWare Project File |
|-|template.py|Template for getting started on a problem|
