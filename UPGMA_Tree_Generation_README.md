# UPGMA Tree Generation - Emma Durrenberger

# Description

Given a fasta file containing three or more DNA sequences, a progressive multiple sequence alignment is performed and used to construct a Newick tree which is outputted as a file named "UPGMA_tree.txt".

# Required Software

C++ GCC 11.4.0

# Packages

Packages used: 
* string
* unistd.h
* iostream
* fstream
* vector
* bits/stdc++.h

# Usage 

Once the program is running, the user will be prompted to enter the fasta file name and the match, mismatch, opening gap and extending gap scores/penalties. 

```bash
g++ UPGMA_Tree_Generation.cpp	
./a.out
```