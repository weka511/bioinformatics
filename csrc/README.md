# C++

| # | Location | Header |Description |
 ---- | -------------------- |  ---------------------  |------------------------------------|
-|Makefile||
DNA|dna.cpp|dna.hpp|[Counting DNA Nucleotides](http://rosalind.info/problems/dna/) |
-|factory.cpp|factory.hpp|This class instantiates the solver for a Rosalind problem 
-|file-adapter.cpp|file-adapter.hpp|Classes that allow Roslind problems to interact with files
-|newick.cpp|newick.hpp|This class parses a string to a Newick tree
-|node.cpp|node.hpp|This class represents a node in a tree
-|problem.cpp|problem.hpp|Abstract classes that serve as parents for a Rosalind problem, and also support input/output
QRTD|qrtd.cpp|qrtd.hpp|[Quartet Distance](http://rosalind.info/problems/qrtd/) [WIP](https://github.com/weka511/bioinformatics/issues/46)|
RNA|rna.cpp|rna.hpp|[Transcribing DNA into RNA](http://rosalind.info/problems/rna/) 
-|rosalind.cpp|rosalind.hpp|Main program, written to solve rosalind problems.
-|test-adapter.cpp|test-adapter.hpp|Classes that allow Rosalind problems to be accessed from unit tests
DNA|test-dna.cpp||Unit tests for DNA
-|test-newick.cpp||Unit tests for newick parser
QRTD|test-qrtd.cpp||Unit tests for QRTD
RNA|test-rna.cpp||Unit tests for rna
-|test-tokenizer.cpp||Unit tests for tokenizer
-|tests.cpp|catch.hpp|[Catch2](https://github.com/catchorg/Catch2) test framework
-|token.cpp|token.hpp|This class is used when parsing a string (e.g. by Newick).  Each Token represents one lexical elemant from string
-|tokenizer.cpp|tokenizer.hpp|This class converts a string to a vector of tokens
