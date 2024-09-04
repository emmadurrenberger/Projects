BLAST Algorithm - Emma Durrenberger

# Description

BLAST (basic local alignment search tool) is a commonly used bioinformatics tool that compares a query sequence to sequences from a database and identifies similar segments between them. Given fasta files containing the database sequences and the query sequence, the sequences are preprocessed by finding all of the kmers (words of length k). Kmers that are found in both the query sequence and a database sequence are used to seed a local alignment. The alignment is extended in both directions of the kmer in the sequences until the highest-scoring alignment for that region is found. The aligned regions, their consensus sequence, and the starting positions in the sequences where the aligned regions occur are outputted. * DNA only *

# Required Software

Python 3.6 or higher

# Modules

Built-in modules used:
* sys
* re

# Packages

Packages used: 
* numpy

All required packages can be installed using the package manager with the following command:

```bash
pip install numpy 
```

# Usage 

The kmer length is set at 28, but can be changed by modifying line 40. Likewise, the match, mismatch, opening gap, and extending gap scores can be changed by modifying lines 181-184.

```bash
python BLAST_Algorithm.py <databases_file> <query_sequence_file>
```