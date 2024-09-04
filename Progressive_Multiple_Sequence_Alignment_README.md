# Progressive Multiple Sequence Alignment - Emma Durrenberger

# Description

Progressive Multiple Sequence Alignment is a multiple sequence alignment technique in which, given three or more DNA sequences, the two most similar sequences are aligned to each other first and a consensus sequence is generated, and then the sequence most similar to the consensus sequence is aligned to it. The process continues until all of the sequences have been aligned.

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

Once the program is running, the user will be prompted to enter the name of the fasta file, and the match, mismatch, opening gap and extending gap scores/penalties. 

```bash
g++ Progressive_Multiple_Sequence_Alignment.cpp	
./a.out
```