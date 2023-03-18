# Bioinformatics
Code inspired by [Bioinformatics Algorithms: an Active Learning Approach](http://bioinformaticsalgorithms.com/) and
by [Rosalind](http://rosalind.info).

NB: functions generally use zero based indexing; [Rosalind](http://rosalind.info) uses 1-based.

I have written this code in Python 3, generally with the most recent release at the time each module was written-- Python 3.11.1 (tags/v3.11.1:a7a450f, Dec  6 2022, 19:58:39) at the time this file was updated. I have used newly added features whenever they might make me more productive. I've been writing code since 1968, during which time computer have got faster, but my brain hasn't. You are welcome to use this code if it is useful to you, subject to the terms of the Licence, but it may need to be adapted if you use an earlier version of Python.

## Bioinformatics Algorithms TextBook Track

### 1. Where in the Genome does DNA replication begin?

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|BA1A|replication.py|[Compute the Number of Times a Pattern Appears in a Text](http://rosalind.info/problems/ba1a/) |
|BA1B|replication.py|[Find the Most Frequent Words in a String](http://rosalind.info/problems/ba1b/) |
|BA1C|replication.py|[Find the Reverse Complement of a String](http://rosalind.info/problems/ba1c/) |
|BA1D|replication.py|[Find All Occurrences of a Pattern in a String](http://rosalind.info/problems/ba1d/) |
|BA1E|replication.py|[Find Patterns Forming Clumps in a String](http://rosalind.info/problems/ba1e/)|
|BA1F|replication.py|[Find a Position in a Genome Minimizing the Skew](http://rosalind.info/problems/ba1f/)|
|BA1G|replication.py|[Compute the Hamming Distance Between Two Strings](http://rosalind.info/problems/ba1g/)|
|BA1H|replication.py|[Find All Approximate Occurrences of a Pattern in a String](http://rosalind.info/problems/ba1h/) |
|BA1I|replication.py|[Find the Most Frequent Words with Mismatches in a String](http://rosalind.info/problems/ba1i/) |
|BA1J|replication.py|[Find Frequent Words with Mismatches and Reverse Complements](http://rosalind.info/problems/ba1j/) |
|BA1K|replication.py|[Generate the Frequency Array of a Strings](http://rosalind.info/problems/ba1f/) |
|BA1L|replication.py|[Implement PatternToNumber](http://rosalind.info/problems/ba1l/) |
|BA1M|replication.py|[Implement NumberToPattern](http://rosalind.info/problems/ba1m/) |
|BA1N|replication.py|[Generate the d-Neighborhood of a String](http://rosalind.info/problems/ba1n/) |
|-|count_kmers.py| Challenge Problems from [Bioinformatics Algorithms](http://bioinformaticsalgorithms.com/): determine number of times each kmer appears in DNA, allowing for reverse complements and d-neighbourhoods.|

### 2. Which DNA Patterns Play the Role of Molecular Clocks?

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|BA2A|sequence.py|[Implement MotifEnumeration](http://rosalind.info/problems/ba2a/) |
|BA2B |sequence.py|[Find a Median String ](http://rosalind.info/problems/ba1j/)
|BA2C |sequence.py|[Find a Profile-most Probable k-mer in a String	](http://rosalind.info/problems/ba2c/)
|BA2D |sequence.py|[Implement GreedyMotifSearch](http://rosalind.info/problems/ba2d/)
|BA2E |sequence.py|[Implement GreedyMotifSearch with Pseudocounts](http://rosalind.info/problems/ba2e/)
|BA2F |sequence.py|[Implement RandomizedMotifSearch](http://rosalind.info/problems/ba2f/)
|BA2G |sequence.py|[Implement GibbsSampler](http://rosalind.info/problems/ba2g/)
|BA2H |sequence.py|[Implement DistanceBetweenPatternAndStrings](http://rosalind.info/problems/ba2h/)

### 3. How Do We Assemble Genomes?  See also  [Genome Assembly](http://rosalind.info/problems/topics/genome-assembly/)

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|BA3A|assemble.py|[Generate the k-mer Composition of a String](http://rosalind.info/problems/ba3a/) |
|BA3B|assemble.py|[Reconstruct a String from its Genome Path](http://rosalind.info/problems/ba3a/)
|BA3C|assemble.py|[Construct the Overlap Graph of a Collection of k-mers](http://rosalind.info/problems/ba3a/)
|BA3D |assemble.py|[Construct the De Bruijn Graph of a String ](http://rosalind.info/problems/ba3a/)
|BA3E |assemble.py|[Construct the De Bruijn Graph of a Collection of k-mers](http://rosalind.info/problems/ba3a/)
|BA3F |assemble.py|[Find an Eulerian Cycle in a Graph](http://rosalind.info/problems/ba3a/)
|BA3G |assemble.py|[Find an Eulerian Path in a Graph](http://rosalind.info/problems/ba3a/)
|BA3H |assemble.py|[Reconstruct a String from its k-mer Composition](http://rosalind.info/problems/ba3a/)
|BA3I |assemble.py|[Find a k-Universal Circular String](http://rosalind.info/problems/ba3a/)
|BA3J |assemble.py|[Reconstruct a String from its Paired Composition](http://rosalind.info/problems/ba3a/)
|asmq|ASMQ.py|[Assessing Assembly Quality with N50 and N75](http://rosalind.info/problems/asmq/) |
|corr|CORR.py|[Error Correction in Reads](http://rosalind.info/problems/corr/) |
|gasm|GASM.py|[Genome Assembly Using Reads](http://rosalind.info/problems/gasm/) |
|grep|GREP.py|[Genome Assembly with Perfect Coverage and Repeats](http://rosalind.info/problems/grep/)|
|KMER|assemble.py|[Generalizing GC-Content](https://rosalind.info/problems/kmer/)|
|long|LONG.py|[Genome Assembly as Shortest Superstring](http://rosalind.info/problems/long/) |
|pcov|PCOV.py|[Genome Assembly with Perfect Coverage](http://rosalind.info/problems/pcov/) |


### 4. How Do We Sequence Antibiotics? See also [Computational Mass Spectrometry](http://rosalind.info/problems/topics/computational-mass-spectrometry/)

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|BA4A|spectrum.py|[Translate an RNA String into an Amino Acid String](http://rosalind.info/problems/ba4a/) |
|BA4B|spectrum.py|[Find Substrings of a Genome Encoding a Given Amino Acid String](http://rosalind.info/problems/ba4b/) |
|BA4C|spectrum.py|[Generate the Theoretical Spectrum of a Cyclic Peptide](http://rosalind.info/problems/ba4c/) |
|BA4D|spectrum.py|[Compute the Number of Peptides of Given Total Mass](http://rosalind.info/problems/ba4d/) |
|BA4E|spectrum.py|[Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum](http://rosalind.info/problems/ba4e/) |
|BA4F|spectrum.py|[Compute the Score of a Cyclic Peptide Against a Spectrum](http://rosalind.info/problems/ba4f/) |
|BA4G|spectrum.py|[Implement LeaderboardCyclopeptideSequencing ](http://rosalind.info/problems/ba4g/) |
|BA4H|spectrum.py|[Generate the Convolution of a Spectrum ](http://rosalind.info/problems/ba4h/) |
|BA4I|spectrum.py|[Implement ConvolutionCyclopeptideSequencing](http://rosalind.info/problems/ba4i/) |
|BA4J|spectrum.py|[Generate the Theoretical Spectrum of a Linear Peptide](http://rosalind.info/problems/ba4j/) |
|BA4K|spectrum.py|[Compute the Score of a Linear Peptide ](http://rosalind.info/problems/ba4k/) |
|BA4L|spectrum.py|[Trim a Peptide Leaderboard](http://rosalind.info/problems/ba4l/) |
|BA4M|ba4m.py spectrum.py| [Solve the Turnpike Problem](http://rosalind.info/problems/ba4m/) |
|-|tyrocidine.py| Challenge Problems from [Bioinformatics Algorithms](http://bioinformaticsalgorithms.com/): sequence Tyrocidine from supplied spectrum (WIP)|
|conv|conv.py+spectrum.py|[Comparing Spectra with the Spectral Convolution](http://rosalind.info/problems/conv/)|
|full|full.py+spectrum.py|[Inferring Peptide from Full Spectrum](http://rosalind.info/problems/full/)|
|prtm|spectrum.py|[Calculating Protein Mass](http://rosalind.info/problems/prtm/) - get_weight(...)|
|prsm|prsm.py+spectrum.py|[Matching a Spectrum to a Protein](http://rosalind.info/problems/prsm/)|
|spec|spec.py+spectrum.py|[Inferring Protein from Spectrum](http://rosalind.info/problems/spec/)|
|sgra|sgra.py+spectrum.py|[Using the Spectrum Graph to Infer Peptides](http://rosalind.info/problems/sgra/)|

### 5. How do we compare Biological Sequences?  See also [Alignment](http://rosalind.info/problems/topics/alignment/)

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
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
|BA5L|BA5L.py align.py|[Align Two Strings Using Linear Space](http://rosalind.info/problems/ba5l/) [WIP](https://github.com/weka511/bioinformatics/issues/29)|
|BA5M|BA5M.py align.py|[Find a Highest-Scoring Multiple Sequence Alignment](http://rosalind.info/problems/ba5m/) |
|BA5N |align.py align.py|[Find a Topological Ordering of a DAG](http://rosalind.info/problems/ba5n/) |
|ctea|ctea.py align.py|[Counting Optimal Alignments](http://rosalind.info/problems/ctea/)|
|edit|edit.py align.py|[Edit Distance](http://rosalind.info/problems/edit/) |
|edta|edta.py align.py|[Edit Distance Alignment](http://rosalind.info/problems/edta/)|
|gaff|gaff.py|[Global Alignment with Scoring Matrix and Affine Gap Penalty](http://rosalind.info/problems/gaff/)|
|gcon|gcon.py|[Global Alignment with Constant Gap Penalty](http://rosalind.info/problems/gcon/)|
|glob|GLOB.py|[Global Alignment with Scoring Matrix](http://rosalind.info/problems/glob/) -- test harness only|
|hamm|rosalind.py|[Counting Point Mutations](http://rosalind.info/problems/hamm/) |
|laff|laff.py |[Local Alignment with Affine Gap Penalty](http://rosalind.info/problems/laff/)|
|loca|loca.py|[Local Alignment with Scoring Matrix ](http://rosalind.info/problems/loca/) |
|mgap|mmgap.py|[Maximizing the Gap Symbols of an Optimal Alignment](http://rosalind.info/problems/mgap/)|
|mult|mult.py align.py|[Multiple Alignment](http://rosalind.info/problems/mult/)|
|oap|oap.py|[Overlap Alignment](http://rosalind.info/problems/oap/) |
|osym|osym.py|[Isolating Symbols in Alignments](http://rosalind.info/problems/osym/) [WIP](https://github.com/weka511/bioinformatics/issues/97)|
|pdst|rosalind.py|[Creating a Distance Matrix ](http://rosalind.info/problems/pdst/) |
|sims|sims.py align.py|[Finding Mutated Motifs](http://rosalind.info/problems/sims/)|
|smgb|smgb.py |[Semiglobal Alignment](http://rosalind.info/problems/smgb/) |
|tran|rosalind.py|[Transitions and Transversions](http://rosalind.info/problems/tran/) |


### 6. Are there fragile regions in the human genome?

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|BA6A|BA6A.py fragile.py|[Implement GreedySorting to Sort a Permutation by Reversals ](http://rosalind.info/problems/ba6a/)|
|BA6B|BA6B.py fragile.py|[Compute the Number of Breakpoints in a Permutation](http://rosalind.info/problems/ba6b/)|
|BA6C|BA6C.py fragile.py|[Compute the 2-Break Distance Between a Pair of Genomes](http://rosalind.info/problems/ba6c/)|
|BA6D|BA6D.py|[Find a Shortest Transformation of One Genome into Another by 2-Breaks](http://rosalind.info/problems/ba6d/) [WIP](https://github.com/weka511/bioinformatics/issues/13)|
|BA6E|fragile.py| [Find All Shared k-mers of a Pair of Strings ](http://rosalind.info/problems/ba6e/)|
|BA6F|BA6F.py fragile.py|[Implement Chromosome to Cycle](http://rosalind.info/problems/ba6f/)|
|BA6G|BA6G.py fragile.py|[Implement Cycle to Chromosome](http://rosalind.info/problems/ba6g/)|
|BA6H|BA6H.py fragile.py|[Implement Coloured Edges](http://rosalind.info/problems/ba6h/)|
|BA6I|BA6I.py fragile.py|[Implement GraphToGenome ](http://rosalind.info/problems/ba6i/)|
|BA6J|BA6J.py fragile.py|[Implement 2-BreakOnGenomeGraph ](http://rosalind.info/problems/ba6j/)|
|BA6K|BA6K.py fragile.py|[Implement 2-BreakOnGenome ](http://rosalind.info/problems/ba6k/)|

### 7. Which animal gave us SARS? See also [Phylogeny](http://rosalind.info/problems/topics/phylogeny/)

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|BA7A|BA7A.py|[Compute Distances Between Leaves](http://rosalind.info/problems/ba7a/)|
|BA7B|BA7B.py|[Limb Length Problem](http://rosalind.info/problems/ba7b/)|
|BA7C|BA7C.py|[Implement Additive Phylogeny](http://rosalind.info/problems/ba7c/) [See also](https://github.com/weka511/bioinformatics/issues/19)|
|BA7D|BA7D.py|[Implement UPGMA](http://rosalind.info/problems/ba7d/)|
|BA7E|BA7E.py|[Implement the Neighbor Joining Algorithm](http://rosalind.info/problems/ba7e/)|
|BA7F|BA7F.py|[Implement SmallParsimony](http://rosalind.info/problems/ba7f/)|
|BA7G|BA7G.py|[Adapt SmallParsimony to Unrooted Trees](http://rosalind.info/problems/ba7g/)|
|alph|alph.py phylogeny.py|[Alignment-Based Phylogeny](http://rosalind.info/problems/alph/)|
|chbp|chbp.py phylogeny.py|[Character-Based Phylogeny](http://rosalind.info/problems/chbp/)|
|cntq|cntq.py phylogeny.py|[Counting Quartets](http://rosalind.info/problems/cntq/)|
|cset|cset.py  phylogeny.py|[Fixing an Inconsistent Character Set](http://rosalind.info/problems/cset/) |
|cstr|cstr.py phylogeny.py|[Creating a Character Table from Genetic Strings](http://rosalind.info/problems/cstr/)|
|ctbl|ctbl.py phylogeny.py|[Creating a Character Table](http://rosalind.info/problems/ctbl/) |
|cunr|phylogeny.py|[Counting Unrooted Binary Trees](http://rosalind.info/problems/cunr/) |
|eubt|eubt.py phylogeny.py|[Enumerating Unrooted Binary Trees](http://rosalind.info/problems/eubt/)|
|inod|-|[Counting Phylogenetic Ancestors](http://rosalind.info/problems/inod/) No code needed|
|mend|mend.py phylogeny.py|[Inferring a Genotype from a Pedigree](http://rosalind.info/problems/mend/) |
|nkew|nwck.py|[Weighting the Tree](http://rosalind.info/problems/nkew/) |
|nwck|nwck.py|[Distances in Trees](http://rosalind.info/problems/nwck/) |
|pdst|rosalind.py|[Creating a Distance Matrix ](http://rosalind.info/problems/pdst/) |
|qrt|qrt.py phylogeny.py|[Quartets](http://rosalind.info/problems/qrt/)|
|qrtd|qrtd.py|[Quartet Distance](http://rosalind.info/problems/qrtd/) [WIP](https://github.com/weka511/bioinformatics/issues/46)|
|root|phylogeny.py|[Counting Rooted Binary Trees ](http://rosalind.info/problems/root/) |
|rsub|rsub.py phylogeny.py|[Identifying Reversing Substitutions](http://rosalind.info/problems/rsub/)|
|sptd|sptb.py phylogeny.py|[Phylogeny Comparison with Split Distance](http://rosalind.info/problems/sptd/) |
|tree|tree.py phylogeny.py|[Completing a Tree ](http://rosalind.info/problems/tree/) |

### 8. How did Yeast become a Winemaker?

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|BA8A|BA8A.py | [Implement FarthestFirstTraversal](http://rosalind.info/problems/ba8a/) |
|BA8B|BA8B.py | [Compute the Squared Error Distortion](http://rosalind.info/problems/ba8b/) |
|BA8C|BA8C.py | [Implement the Lloyd Algorithm for k-Means Clustering](http://rosalind.info/problems/ba8c/) |
|BA8D|BA8D.py | [Implement the Soft k-Means Clustering Algorithm](http://rosalind.info/problems/ba8d/) |
|BA8E|BA8E.py | [Implement Hierarchical Clustering](http://rosalind.info/problems/ba8e/) |

### 9. How do we locate disease causing mutations?

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
||snp.py|Code for: How do we locate disease causing mutations?|
|BA9A|BA9A.py snp.py| [Construct a Trie from a Collection of Patterns](http://rosalind.info/problems/ba9a/) |
|BA9B|BA9B.py snp.py| [Implement TrieMatching](http://rosalind.info/problems/ba9b/) |
|BA9C|BA9C.py snp.py| [Construct the Suffix Tree of a String](http://rosalind.info/problems/ba9c/)|
|BA9D|BA9D.py| [Find the Longest Repeat in a String](http://rosalind.info/problems/ba9d/)|
|BA9E|BA9E.py snp.py| [Find the Longest Substring Shared by Two Strings](http://rosalind.info/problems/ba9e/)|
|BA9F|BA9F.py| [Find the Shortest Non-Shared Substring of Two Strings](http://rosalind.info/problems/ba9f/) [WIP](https://github.com/weka511/bioinformatics/issues/73)|
|BA9G|BA9G.py snp.py| [Construct the Suffix Array of a String](http://rosalind.info/problems/ba9g/) |
|BA9H|BA9H.py| [Pattern Matching with the Suffix Array ](http://rosalind.info/problems/ba9h/) |
|BA9I|BA9I.py snp.py| [Construct the Burrows-Wheeler Transform of a String](http://rosalind.info/problems/ba9i/) |
|BA9J|BA9J.py snp.py| [Reconstruct a String from its Burrows-Wheeler Transform](http://rosalind.info/problems/ba9j/) |
|BA9K|BA9K.py snp.py| [Generate the Last-to-First Mapping of a String](http://rosalind.info/problems/ba9k/) |
|BA9L|BA9L.py snp.py| [Implement BWMatching ](http://rosalind.info/problems/ba9l/)|
|BA9M|BA9M.py snp.py| [Implement BetterBWMatching](http://rosalind.info/problems/ba9m/)|
|BA9N|BA9N.py snp.py| [Find All Occurrences of a Collection of Patterns in a String ](http://rosalind.info/problems/ba9n/) |
|BA9O|BA9O.py| [Find All Approximate Occurrences of a Collection of Patterns in a String](http://rosalind.info/problems/ba9o/) [WIP](https://github.com/weka511/bioinformatics/issues/79)|
|BA9P|BA9P.py| [Implement Tree Colouring](http://rosalind.info/problems/ba9p/)|
|BA9Q|BA9Q.py snp.py| [Construct the Partial Suffix Array of a String](http://rosalind.info/problems/ba9q/)|
|BA9R|BA9R.py| [Construct a Suffix Tree from a Suffix Array](http://rosalind.info/problems/ba9r/) [WIP](https://github.com/weka511/bioinformatics/issues/81)|

### 10. Why have biologists still not developed an HIV vaccine?

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|BA10A|BA10A.py hmm.py| [Compute the Probability of a Hidden Path](http://rosalind.info/problems/ba10a/) |
|BA10B|BA10B.py hmm.py| [Compute the Probability of an Outcome Given a Hidden Path](http://rosalind.info/problems/ba10b/) |
|BA10C|BA10C.py hmm.py| [Implement the Viterbi Algorithm](http://rosalind.info/problems/ba10c/) |
|BA10D|BA10D.py hmm.py| [Compute the Probability of a String Emitted by an HMM ](http://rosalind.info/problems/ba10d/)|
|BA10E|BA10E.py hmm.py| [Construct a Profile HMM](http://rosalind.info/problems/ba10e/) |
|BA10F|BA10F.py hmm.py| [Construct a Profile HMM with Pseudocounts](http://rosalind.info/problems/ba10f/)|
|BA10G|BA10G.py hmm.py| [Perform a Multiple Sequence Alignment with a Profile HMM ](http://rosalind.info/problems/ba10g/) [WIP](https://github.com/weka511/bioinformatics/issues/85)|
|BA10H|BA10H.py hmm.py| [Estimate the Parameters of an HMM ](http://rosalind.info/problems/ba10h/) |
|BA10I|BA10I.py hmm.py| [Implement Viterbi Learning ](http://rosalind.info/problems/ba10i/) |
|BA10J|BA10J.py hmm.py| [Solve the Soft Decoding Problem ](http://rosalind.info/problems/ba10j/) |
|BA10K|BA10K.py hmm.py| [Implement Baum-Welch Learning](http://rosalind.info/problems/ba10k/)|

### 11.  Was T.rex just a big chicken?

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|BA11A|BA11A.py spectrum.py| [Spectrum Graph Construction](http://rosalind.info/problems/ba11a/)|
|BA11B|BA11B.py spectrum.py| [Implement DecodingIdealSpectrum](http://rosalind.info/problems/ba11b/) |
|BA11C|BA11C.py spectrum.py| [Convert a Peptide into a Peptide Vector ](http://rosalind.info/problems/ba11c/) |
|BA11D|BA11D.py spectrum.py| [Convert a Peptide Vector into a Peptide ](http://rosalind.info/problems/ba11d/) |
|BA11E|BA11E.py spectrum.py| [Sequence a Peptide](http://rosalind.info/problems/ba11e/) [WIP](https://github.com/weka511/bioinformatics/issues/91)|
|BA11F|BA11F.py spectrum.py| [Find a Highest-Scoring Peptide in a Proteome against a Spectrum](http://rosalind.info/problems/ba11f/) [WIP](https://github.com/weka511/bioinformatics/issues/92)|
|BA11G|BA11G.py spectrum.py| [Implement PSMSearch](http://rosalind.info/problems/ba11g/) [WIP](https://github.com/weka511/bioinformatics/issues/93)|
|BA11H|BA11H.py spectrum.py| [Compute the Size of a Spectral Dictionary](http://rosalind.info/problems/ba11h/) [WIP](https://github.com/weka511/bioinformatics/issues/94)|
|BA11I|BA11I.py spectrum.py| [Compute the Probability of a Spectral Dictionary](http://rosalind.info/problems/ba11i/) [WIP](https://github.com/weka511/bioinformatics/issues/95)|
|BA11J|BA11J.py spectrum.py| [Find a Highest-Scoring Modified Peptide against a Spectrum](http://rosalind.info/problems/ba11j/) [WIP](https://github.com/weka511/bioinformatics/issues/96)|



## [Combinatorics](http://rosalind.info/problems/topics/combinatorics/)

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
||combinatorics.py|[Combinatorics](http://rosalind.info/problems/topics/combinatorics/)|
|aspc|rosalind.py|[Introduction to Alternative Splicing](http://rosalind.info/problems/aspc/) |
|cat|cat.py combinatorics.py|[Catalan Numbers and RNA Secondary Structures ](http://rosalind.info/problems/cat/)|
|cntq|cntq.py |[Counting Quartets](http://rosalind.info/problems/cntq/)|
|cunr|phylogeny.py|[Counting Unrooted Binary Trees](http://rosalind.info/problems/cunr/) |
|eubt|eubt.py phylogeny.py|[Enumerating Unrooted Binary Trees](http://rosalind.info/problems/eubt/)|
|fibd|rosalind.py|[Mortal Fibonacci Rabbits](http://rosalind.info/problems/fibd/) |
|fib|rosalind.py|[Rabbits and Recurrence Relations](http://rosalind.info/problems/fib/) |
|motz|motz.py combinatorics.py|[Motzkin Numbers and RNA Secondary Structures](http://rosalind.info/problems/motz/) |
|mrep|mrep.py|[Identifying Maximal Repeats ](http://rosalind.info/problems/mrep/)|
|mrna|rosalind.py|[Inferring mRNA from Protein](http://rosalind.info/problems/mrna/) |
|orf|orf.py|[Open Reading Frames](http://rosalind.info/problems/orf/)|
|pmch|pmch.py|[Perfect Matchings and RNA Secondary Structures](http://rosalind.info/problems/pmch/) |
|pper|rosalind.py|[Partial Permutations](http://rosalind.info/problems/pper/) |
|rnas|rnas.py|[Wobble Bonding and RNA Secondary Structures](http://rosalind.info/problems/rnas/) [WIP](https://github.com/weka511/bioinformatics/issues/70)|
|root|phylogeny.py|[Counting Rooted Binary Trees ](http://rosalind.info/problems/root/) |



## [Genome Rearrangements](http://rosalind.info/problems/topics/genome-rearrangements/)

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|lgis|lgis.py|[Longest Increasing Subsequence](http://rosalind.info/problems/lgis/) |
|rear|rear.py fragile.py|[Reversal distance](http://rosalind.info/problems/rear/) |
|sort|sort.py fragile.py|[Sorting by Reversals](http://rosalind.info/problems/sort/) |

## [Graphs](http://rosalind.info/problems/topics/graphs/) & [Graph Algorithms](http://rosalind.info/problems/topics/graph-algorithms/)

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|2sat|sat.py graphs.py|[2-Satisfiability](http://rosalind.info/problems/2sat/)|
|bf|bf.py graphs.py|[Bellman-Ford Algorithm](http://rosalind.info/problems/bf/)|
|bfs|bfs.py graphs.py|[Breadth First search](http://rosalind.info/problems/bfs/) |
|bins|BINS.py|[Binary Search](http://rosalind.info/problems/bins/) |
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
|lrep|lrep.py|[Finding the Longest Multiple Repeat](http://rosalind.info/problems/lrep/) [WIP](https://github.com/weka511/bioinformatics/issues/51)|
|nwc|nwc.py graphs.py|[Negative Weight Cycle](http://rosalind.info/problems/nwc/)|
|scc|scc.py graphs.py|[Strongly Connected Components](http://rosalind.info/problems/scc/)|
|sc|sc.py graphs.py|[Semi-Connected Graph](http://rosalind.info/problems/sc/)|
|sdag|sdag.py graphs.py|[Shortest Paths in DAG](http://rosalind.info/problems/sdag/)|
|sq|sq.py graphs.py|[Square in a Graph](http://rosalind.info/problems/sq/)|
|suff|suff.py  graphs.py|[Creating a Suffix Tree](http://rosalind.info/problems/suff/) |
|tree|tree.py phylogeny.py|[Completing a Tree ](http://rosalind.info/problems/tree/) |
|trie|trie.py snp.py|[Introduction to Pattern Matching](http://rosalind.info/problems/trie/) |
|ts|ts.py align.py|[Topological sort](http://rosalind.info/problems/ts/)|

## [Heredity](http://rosalind.info/problems/topics/heredity/)
| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|iev|rosalind.py| [Calculating Expected Offspring](http://rosalind.info/problems/iev/)|
|indq|rosalind.py|[Independent Segregation of Chromosomes](http://rosalind.info/problems/indq/) |
|iprb|rosalind.py|[Mendel's First Law](http://rosalind.info/problems/iprb/) |
|lia|rosalind.py|[Independent Alleles](http://rosalind.info/problems/lia/) |
|mend|mend.py phylogeny.py|[Inferring a Genotype from a Pedigree](http://rosalind.info/problems/mend/) |
|sexl|rosalind.py|[Sex-Linked Inheritance ](http://rosalind.info/problems/sexl/) |



## [Population Dynamics](http://rosalind.info/problems/topics/population-dynamics/)

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|afrq|rosalind.py|[Counting Disease Carriers](http://rosalind.info/problems/afrq/) |
|foun|rosalind.py|[The Founder Effect and Genetic Drift](http://rosalind.info/problems/foun/) |
|wfmd|rosalind.py|[The Wright-Fisher Model of Genetic Drift](http://rosalind.info/problems/wfmd/) |

## [Probability](http://rosalind.info/problems/topics/probability/)

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|ebin|rosalind.py|[Wright-Fisher's Expected Behavior](http://rosalind.info/problems/ebin/) |
|eval|rosalind.py|[Expected Number of Restriction Sites](http://rosalind.info/problems/eval/) |
|mend|mend.py phylogeny.py|[Inferring a Genotype from a Pedigree](http://rosalind.info/problems/mend/) |
|rstr|rosalind.py|[Matching Random Motifs](http://rosalind.info/problems/rstr/) |

## [Sequence Analysis](http://rosalind.info/problems/topics/sequence-analysis/)
| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|frmt|FRMT.py|[Data Formats](http://rosalind.info/problems/frmt/)|
|gbk|GBK.py|[GenBank Introduction](http://rosalind.info/problems/gbk/) |
|orfr|orfr.py|[Finding Genes with Open Reading Frames](http://rosalind.info/problems/orfr/) |
|ptra|PTRA.py|[Protein Translation](http://rosalind.info/problems/ptra/) |

## [Sorting](http://rosalind.info/problems/topics/sorting/)

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
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

## Strings & [String Algorithms](http://rosalind.info/problems/topics/string-algorithms/)
| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|cons|rosalind.py|[Consensus and Profile](http://rosalind.info/problems/cons/) |
|dna|dna.py rosalind.py|[Counting DNA Nucleotides](http://rosalind.info/problems/dna/) |
|gc|rosalind.py|[Computing GC Content](http://rosalind.info/problems/gc/) |
|itwv|itwv.py align.py|[Finding Disjoint Motifs in a Gene ](http://rosalind.info/problems/itwv/)|
|kmer|rosalind.py|[Generalizing GC-Content](http://rosalind.info/problems/kmer/) |
|kmp|kmp.py|[Speeding Up Motif Finding](http://rosalind.info/problems/kmp/) |
|ksim|ksim.py|[Finding All Similar Motifs](http://rosalind.info/problems/ksim/) [WIP](https://github.com/weka511/bioinformatics/issues/62)|
|lcsm|rosalind.py|[Finding a Shared Motif](http://rosalind.info/problems/lcsm/) |
|lcsq|lcsq.py|[Finding a Shared Spliced Motif](http://rosalind.info/problems/lcsq/) |
|ling|ling.py|[Linguistic Complexity of a Genome](http://rosalind.info/problems/ling/)|
||ukkonen.py|[µTuX's implementation of Ukkonen's algorithm](https://github.com/mutux/Ukkonen-s-Suffix-Tree-Algorithm/blob/master/suffixtree.py)|
|mmch|mmch.py|[Maximum Matchings and RNA Secondary Structures](http://rosalind.info/problems/mmch/) |
|prot|spectrum.py|[Translating RNA into Protein](http://rosalind.info/problems/prot/) |
|revc|rosalind.py| [Complementing a Strand of DNA](http://rosalind.info/problems/revc/)|
|revp|revp.py rosalind.py| [Locating Restriction Sites](http://rosalind.info/problems/revp/)|
|rna|rna.py rosalind.py|[Transcribing DNA into RNA](http://rosalind.info/problems/rna/) |
|scsp|scsp.py|[Interleaving Two Motifs](http://rosalind.info/problems/scsp/)|
|splc|spectrum.py|[RNA Splicing](http://rosalind.info/problems/splc/) |
|sseq|rosalind.py|[Finding a Spliced Motif](http://rosalind.info/problems/sseq/) |
|subs|rosalind.py|[Finding a Motif in DNA](http://rosalind.info/problems/subs/) |



## [Set Theory](http://rosalind.info/problems/topics/set-theory/)
| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|pdpl|pdpl|[Creating a Restriction Map](http://rosalind.info/problems/pfpl/)|
|seto|-|[Set operations](http://rosalind.info/problems/seto/) Solved using Python shell|
|sset|-|[Counting Subsets](http://rosalind.info/problems/sset/) Solved with a single line in Python shell: 2**n%1000000|



## Other Problems

| # | Location | Description |
| ---- | --------------- | ------------------------------------------------|
|fibo|fibo.py|[Fibonacci numbers](http://rosalind.info/problems/fibo/) |



## Support

| Location | Description |
| --------------- | ------------------------------------------------|
|bioinformatics.bib|Bibliography|
|helpers.py|Utilities for formatting output, parsing input, etc |
|LICENSE|License Agreement|
|newick.py|Parser for files in [Newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html) |
|README.md|This file|
|rosalind.py|Shared code|
|rosalind.wpr|WingWare Project File |
|template.py|Template for getting started on a problem|
