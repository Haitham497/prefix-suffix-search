# prefix-suffix-search
A data-structures–based implementation of Prefix Trees (Tries) and Suffix Trees for fast pattern matching in long biological sequences. This project demonstrates how tree-based indexing can efficiently search for short sequences—such as genes, motifs, or regulatory elements—within large genomes.

DNA Subsequence Search Using Prefix and Suffix Trees
This project uses two types of tree data structures — Prefix Trees (Tries) and Suffix Trees — to search for DNA sequences inside a genome. It helps in searching faster and using less memory.
Methods:
Prefix Tree (Trie): Stores short DNA pieces (kmers) from the genome and searches for matches. It allows small errors like substitutions.
Suffix Tree: Stores all possible endings of the genome and searches for exact matches of any sequence.
Files Included:
prefixTree.cpp, prefixTree.h, and main.cpp.
suffixTree.cpp, suffixTree.h, and main.cpp.
covid_partial.fasta: Genome file used for the suffix tree (only a small part of the COVID-19 genome).
How It Works:
Prefix Tree
-	Breaks the genome into small fixed-length pieces (kmers).
-	Stores them in a trie.
-	You can search for sequences with or without small errors.
-	Works well for short reads and fast lookup.
Suffix Tree
-	Builds a tree from every ending of the genome.
-	Used to search for any sequence exactly.
-	Works well for full matches.
•	Current version builds a basic tree without checking for duplicate paths, so it's slower and uses more memory.

How to Run
Open the project in an IDE like Visual Studio or use a C++ compiler.
Make sure the genome file is in the same folder.
Build and run the program.

Notes
The suffix tree uses a partial COVID-19 genome to save time and resources.
The prefix tree is faster and more efficient for short reads.
The suffix tree can be improved by avoiding repeated paths while building.

