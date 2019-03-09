# Bioinformatics

Code inspired by [Finding Hidden Messages in DNA (Bioinformatics I)](https://class.coursera.org/hiddenmessages-003) and 
from [Rosalind](http://rosalind.info). 
NB: functions generally use zero based indexing; Rosalind uses 1-based.

## Textbook track

| # | Location | Description |
| ---- | -------------------------- | ------------------------------------------------|
|1|-| Where in the Genome does DNA replication begin?|
|BA1A|rosalind.py|	[Compute the Number of Times a Pattern Appears in a Text](http://rosalind.info/problems/ba1a/) |
|BA1B|rosalind_old.py|	[Find the Most Frequent Words in a String](http://rosalind.info/problems/ba1b/) |
|BA1C|rosalind_old.py|	[Find the Reverse Complement of a String](http://rosalind.info/problems/ba1c/) |
|BA1D|rosalind_old.py|	[Find All Occurrences of a Pattern in a String](http://rosalind.info/problems/ba1d/) |
|BA1E|rosalind_old.py|	[Find Patterns Forming Clumps in a String](http://rosalind.info/problems/ba1e/) |
|BA1F|rosalind_old.py|	[Find a Position in a Genome Minimizing the Skew](http://rosalind.info/problems/ba1f/) |
|BA1G|rosalind_old.py|	[Compute the Hamming Distance Between Two Strings](http://rosalind.info/problems/ba1g/) |
|BA1H|rosalind_old.py|	[Find All Approximate Occurrences of a Pattern in a String](http://rosalind.info/problems/ba1h/) |
|BA1I|rosalind_old.py|	[Find the Most Frequent Words with Mismatches in a String](http://rosalind.info/problems/ba1i/) |
|BA1J|rosalind_old.py|	[Find Frequent Words with Mismatches and Reverse Complements](http://rosalind.info/problems/ba1j/) |
|BA1K|rosalind_old.py|	[Generate the Frequency Array of a Strings](http://rosalind.info/problems/ba1f/) |
|BA1L|rosalind_old.py|	[Implement PatternToNumber](http://rosalind.info/problems/ba1l/) |
|BA1M|rosalind_old.py|	[Implement NumberToPattern](http://rosalind.info/problems/ba1m/) |
|BA1N|rosalind_old.py|[Generate the d-Neighborhood of a String](http://rosalind.info/problems/ba1n/) |
|2|-|How do we sequence Antibiotics?|
|BA2A|rosalind_old.py| 	[Implement MotifEnumeration](http://rosalind.info/problems/ba2a/) |
|BA2B |rosalind_old.py|	[Find a Median String ](http://rosalind.info/problems/ba1j/)
|BA2C |rosalind_old.py|	[Find a Profile-most Probable k-mer in a String	](http://rosalind.info/problems/ba2c/)
|BA2D |rosalind_old.py|	[Implement GreedyMotifSearch](http://rosalind.info/problems/ba2d/)
|BA2E |rosalind_old.py|[Implement GreedyMotifSearch with Pseudocounts](http://rosalind.info/problems/ba2e/)
|BA2F |rosalind_old.py|	[Implement RandomizedMotifSearch](http://rosalind.info/problems/ba2f/)	
|BA2G |rosalind_old.py|	[Implement GibbsSampler](http://rosalind.info/problems/ba2g/)
|BA2H |rosalind_old.py|	[Implement DistanceBetweenPatternAndStrings](http://rosalind.info/problems/ba2h/) 
|3|-|Which DNA Patterna play the Role of Molecular Clocks?|
|BA3A|rosalind_old.py|	[Generate the k-mer Composition of a String](http://rosalind.info/problems/ba3a/) |
|BA3B|rosalind_old.py|	[Reconstruct a String from its Genome Path](http://rosalind.info/problems/ba3a/)
|BA3C|rosalind_old.py|[Construct the Overlap Graph of a Collection of k-mers](http://rosalind.info/problems/ba3a/)
|BA3D |rosalind_old.py|	[Construct the De Bruijn Graph of a String ](http://rosalind.info/problems/ba3a/)
|BA3E |rosalind_old.py|	[Construct the De Bruijn Graph of a Collection of k-mers](http://rosalind.info/problems/ba3a/)
|BA3F |rosalind_old.py|	[Find an Eulerian Cycle in a Graph](http://rosalind.info/problems/ba3a/)
|BA3G |rosalind_old.py|	[Find an Eulerian Path in a Graph](http://rosalind.info/problems/ba3a/)
|BA3H |rosalind_old.py|	[Reconstruct a String from its k-mer Composition](http://rosalind.info/problems/ba3a/)
|BA3I |rosalind_old.py|	[Find a k-Universal Circular String](http://rosalind.info/problems/ba3a/)
|BA3J |rosalind_old.py|	[Reconstruct a String from its Paired Composition](http://rosalind.info/problems/ba3a/)
|4|-|How do we assemble Genomes?|
|BA4A|rosalind_old.py|[Translate an RNA String into an Amino Acid String](http://rosalind.info/problems/ba4a/) |
|BA4B|rosalind_old.py|[Find Substrings of a Genome Encoding a Given Amino Acid String](http://rosalind.info/problems/ba4b/) |
|BA4C|rosalind_old.py|[Generate the Theoretical Spectrum of a Cyclic Peptide](http://rosalind.info/problems/ba4c/) |
|BA4D|rosalind_old.py|[Compute the Number of Peptides of Given Total Mass](http://rosalind.info/problems/ba4d/) |
|BA4E|rosalind_old.py|	[Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum](http://rosalind.info/problems/ba4e/) |
|BA4F|rosalind_old.py|	[Compute the Score of a Cyclic Peptide Against a Spectrum](http://rosalind.info/problems/ba4f/) |
|BA4G|rosalind_old.py|	[Implement LeaderboardCyclopeptideSequencing ](http://rosalind.info/problems/ba4g/) |
|BA4H|rosalind_old.py|	[Generate the Convolution of a Spectrum ](http://rosalind.info/problems/ba4h/) |
|BA4I|rosalind_old.py|	[Implement ConvolutionCyclopeptideSequencing](http://rosalind.info/problems/ba4i/) |
|BA4J|rosalind_old.py|	[Generate the Theoretical Spectrum of a Linear Peptide](http://rosalind.info/problems/ba4j/) |
|BA4K|rosalind_old.py|	[Compute the Score of a Linear Peptide ](http://rosalind.info/problems/ba4k/) |
|BA4L|rosalind_old.py|	[Trim a Peptide Leaderboard](http://rosalind.info/problems/ba4l/) |
|BA4M|ba4m.py| [Solve the Turnpike Problem](http://rosalind.info/problems/ba4m/) |
|5|-|How do we compare Biological Sequences?|
|BA5A|rosalind_old.py|[Find the Minimum Number of Coins Needed to Make Change](http://rosalind.info/problems/ba5a/)|
|BA5B|rosalind_old.py|[Find the Length of a Longest Path in a Manhattan-like Grid](http://rosalind.info/problems/ba5b/)|
|BA5C|rosalind_old.py|[Find a Longest Common Subsequence of Two Strings](http://rosalind.info/problems/ba5c/)|
|BA5D|rosalind_old.py|[Find the Longest Path in a DAG](http://rosalind.info/problems/ba5d/)|
|BA5E|rosalind_old.py|[Find a Highest-Scoring Alignment of Two Strings](http://rosalind.info/problems/ba5e/)|
|BA5F|rosalind_old.py|[Find a Highest-Scoring Local Alignment of Two Strings](http://rosalind.info/problems/ba5f/)|
|BA5G | edit.py | [Compute the Edit Distance Between Two Strings](http://rosalind.info/problems/ba5g/) |
|BA5H| BA5H.py | [Find a Highest-Scoring Fitting Alignment of Two Strings](http://rosalind.info/problems/ba5h/) |   
|BA5I|BA5I.py|[Find a Highest-Scoring Overlap Alignment of Two Strings](http://rosalind.info/problems/ba5i/)|    
|BA5J|BA5J.py|[Align Two Strings Using Affine Gap Penalties](http://rosalind.info/problems/ba5j/)|    
|BA5K|TBP|[Find a Middle Edge in an Alignment Graph in Linear Space](http://rosalind.info/problems/ba5k/) |
|BA5L|TBP|[Align Two Strings Using Linear Space](http://rosalind.info/problems/ba5l/) |    
|BA5M|TBP|[Find a Highest-Scoring Multiple Sequence Alignment](http://rosalind.info/problems/ba5m/) | 
|BA5N |align.py|[Find a Topological Ordering of a DAG](http://rosalind.info/problems/ba5n/) | 
|6|fragile.py|Are there fragile regions in the human genome?|
|BA6A|BA6A.py| [Implement GreedySorting to Sort a Permutation by Reversals ](http://rosalind.info/problems/ba6a/)|
|BA6B| BA6B.py | [Compute the Number of Breakpoints in a Permutation](http://rosalind.info/problems/ba6b/)|
|BA6C| BA6C.py| [Compute the 2-Break Distance Between a Pair of Genomes](http://rosalind.info/problems/ba6c/)|
|BA6D| TBP| [Find a Shortest Transformation of One Genome into Another by 2-Breaks](http://rosalind.info/problems/ba6d/)|
|BA6E| fragile.py| [Find All Shared k-mers of a Pair of Strings ](http://rosalind.info/problems/ba6e/)|
|BA6F| BA6F.py| [Implement Chromosome to Cycle](http://rosalind.info/problems/ba6f/)|
|BA6G| BA6G.py| [Implement Cycle to Chromosome](http://rosalind.info/problems/ba6g/)|
|BA6H| BA6H.py| [Implement Coloured Edges](http://rosalind.info/problems/ba6h/)|
|BA6I| BA6I.py| [Implement GraphToGenome ](http://rosalind.info/problems/ba6i/)|
|BA6J| BA6J.py| [Implement 2-BreakOnGenomeGraph ](http://rosalind.info/problems/ba6j/)|
|BA6K| BA6K| [Implement 2-BreakOnGenome ](http://rosalind.info/problems/ba6k/) (WIP)|
|7|-|Which animal gave us SARS?|
| BA7A | BA7A.py | [Compute Distances Between Leaves](http://rosalind.info/problems/ba7a/) |
| BA7B |BA7B.py | [Limb Length Problem](http://rosalind.info/problems/ba7b/) |
| BA7C| BA7C.py | [Implement Additive Phylogeny](http://rosalind.info/problems/ba7c/) |
| BA7D|BA7D.py | [Implement UPGMA](http://rosalind.info/problems/ba7d/) |
| BA7E|BA7E.py | [Implement the Neighbor Joining Algorithm](http://rosalind.info/problems/ba7e/) |
| BA7F|BA7F.py | [Implement SmallParsimony](http://rosalind.info/problems/ba7f/) |
| BA7G|BA7.py | [Adapt SmallParsimony to Unrooted Trees](http://rosalind.info/problems/ba7g/) |
|8|-|How did Yeast become a Winemaker?|
|BA8A|BA8A.py | [Implement FarthestFirstTraversal](http://rosalind.info/problems/ba8a/) |
|BA8B|BA8B.py | [Compute the Squared Error Distortion](http://rosalind.info/problems/ba8b/) |
|BA8C|BA8C.py | [Implement the Lloyd Algorithm for k-Means Clustering](http://rosalind.info/problems/ba8c/) |
|BA8D|BA8D.py | [Implement the Soft k-Means Clustering Algorithm](http://rosalind.info/problems/ba8d/) |
|BA8E|BA8E.py | [Implement Hierarchical Clustering](http://rosalind.info/problems/ba8e/) |
|9|-|How do we locate disease causing mutations?|
|BA9A|BA9A.py | [Construct a Trie from a Collection of Patterns](http://rosalind.info/problems/ba9a/) |
|BA9B|BA9B.py | [Implement TrieMatching](http://rosalind.info/problems/ba9b/) |
|10|-|Why have biologists still not developed an HIV vaccine?|
|11|-|Was T.rex just a big chicken?|

## Other Rosalind problems

| Name | Description |
| -------------------------- | ------------------------------------------------|
| 2sum.py | [2SUM](http://rosalind.info/problems/2sum/) |
| 3sum.py | [3SUM](http://rosalind.info/problems/3sum/) |
| ASMQ.py | 	[Assessing Assembly Quality with N50 and N75](http://rosalind.info/problems/asmq/) |
| bfs.py | [Beadth First search](http://rosalind.info/problems/bfs/) |
| bip.py | [Testing Bipartiteness](http://rosalind.info/problems/bip/) |
| BINS.py | [Binary Search](http://rosalind.info/problems/bfs/) |
| cat.pt | [Catalan Numbers and RNA Secondary Structures ](http://rosalind.info/problems/cat/) (WIP)|
| cc.py |	[Connected Components](http://rosalind.info/problems/cc/) |
| CORR.py | [Error Correction in Reads](http://rosalind.info/problems/corr/) | 
| cstr.py | [Creating a Character Table from Genetic Strings](http://rosalind.info/problems/cstr/)|
| ctbl.py | [Creating a Character Table](http://rosalind.info/problems/ctbl/) |
| dag.py | [Testing Acyclicity](http://rosalind.info/problems/dag/) |
| DBRU.py |	[Constructing a De Bruijn Graph](http://rosalind.info/problems/dbru/) |
| DDEG.py | [Double-Degree Array](http://rosalind.info/problems/ddeg/) |
| DEG.py | 	[Degree Array](http://rosalind.info/problems/deg/) |
| edit.py |	[Edit Distance](http://rosalind.info/problems/edit/) |
| edta.py | [Edit Distance Alignment](http://rosalind.info/problems/edta/)|
| FRMT.py | [Data Formats](http://rosalind.info/problems/frmt/)|
| GASM.py | [Genome Assembly Using Reads](http://rosalind.info/problems/gasm/) |
| GBK.py | 	[GenBank Introduction](http://rosalind.info/problems/gbk/) |
| gcon.py | [Global Alignment with Constant Gap Penalty](http://rosalind.info/problems/gcon/)|
| GLOB.py | [Global Alignment with Scoring Matrix](http://rosalind.info/problems/glob/) -- test harness only|
| GREP.py | [Genome Assembly with Perfect Coverage and Repeats](http://rosalind.info/problems/grep/)|
|HEA.py| [Building a Heap](http://rosalind.info/problems/hea/) and [Heap Sort](http://rosalind.info/problems/hs/)
| INS.py | 	[Insertion Sort](http://rosalind.info/problems/ins/) |
| INV.py | 	[Counting Inversions](http://rosalind.info/problems/inv/) |
| kmp.py | [Speeding Up Motif Finding](http://rosalind.info/problems/kmp/) |
| lcsq.py | [Finding a Shared Spliced Motif](http://rosalind.info/problems/lcsq/) |
| ling.py | [Linguistic Complexity of a Genome](http://rosalind.info/problems/ling/) (WIP)|
| loca.py | [Local Alignment with Scoring Matrix ](http://rosalind.info/problems/loca/) |
| LONG.py | [Genome Assembly as Shortest Superstring](http://rosalind.info/problems/long/) |
| maj.py | [Majority Element](http://rosalind.info/problems/maj/) ||
| MED.py | 	[Median](http://rosalind.info/problems/med/)|
| mend.py |  [Inferring a Genotype from a Pedigree](http://rosalind.info/problems/mend/) |
| MER.py | 	[Merge Two Sorted Arrays](http://rosalind.info/problems/mer/) |
| MS.py | 	[Mergesort](http://rosalind.info/problems/ms/)|
| nwck.py | [Distances in Trees](http://rosalind.info/problems/nwck/) |
| oap.py | [Overlap Alignment](http://rosalind.info/problems/oap/) (wip)|
| orfr.py | [Finding Genes with Open Reading Frames](http://rosalind.info/problems/orfr/) |
| PAR.py | 	 [2-Way Partition](http://rosalind.info/problems/par/) |
| PAR3.py |  [3-Way Partition](http://rosalind.info/problems/par3/) |
| PCOV.py |  [Genome Assembly with Perfect Coverage](http://rosalind.info/problems/pcov/) |
| pmch.py| 	[Perfect Matchings and RNA Secondary Structures](http://rosalind.info/problems/pmch/) |
| ps.py | 	[Partial sort](http://rosalind.info/problems/ps/)|
| PTRA.py |  [Protein Translation](http://rosalind.info/problems/ptra/) |
| QS.py | 	[Quicksort](http://rosalind.info/problems/qs/)|
| rear.py | [Reversal distance](http://rosalind.info/problems/rear/) (WIP)|
| ts.py | [Topological sort](http://rosalind.info/problems/ts/) test harness only--see BA5N |


## Support

| Name | |Description |
| -------------------------- | ---- | ------------------------------------------------|
| align.py || Common functions for alignment |
| helpers.py || Utilities for formatting output, parsing input, etc |
| LICENSE ||	License Agreement|
| newick.py || Parser for files in [Newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html) |
| README.md ||This file|
| rosalind.py || Shared code|
|  |revc| [Complementing a Strand of DNA](http://rosalind.info/problems/revc/)|
|  |subs|[Finding a Motif in DNA](http://rosalind.info/problems/subs/) |
|  |trie|[Introduction to Pattern Matching](http://rosalind.info/problems/trie/) |
| rosalind_old.py || Code from the early days of project, including much of the textbook track|
| rosalind.wpr || WingWare Project File |

## Other stuff

| Name | Description |
| ------------------------- | ------------------------------------------------|
| count-kmer-occurrences.py | Determine number of times each kmer appears in DNA, allowing for reverse complements and d-neighbourhoods |
