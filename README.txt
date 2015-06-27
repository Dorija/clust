
To use the Reducer, first create a binary executable file. The file can be created
typing a command:
g++ -std=c++0x -fopenmp src/reducer.cpp -o Reducer

To create a binary executable file for the Searcher, type a command:
g++ -std=c++0x src/searcher.cpp -o Searcher

The output of the Reducer are two files, the file with centroids (reduced database)
and the file with clusters. The Searcher prints the output on the standard output. First
it prints the name of the sequence which we want to identify followed by colon (:) and
the list of all found sequences.

EXAMPLE:
Here is a example of running the Reducer and the Searcher. Please go through an
example before running on you own data.

Fist run the Reducer on test.fasta:
./Reducer --in data/test.fasta --out centroids.fasta --clusters cluster.txt --identity 0.9

To identify a new sequence (seq.fasta):
./Searcher --in data/seq.fasta --dataset centroids.fasta --clusters cluster.txt --identity 0.9

The output should be as follows:
test1:
 seq32
test2:
 seq24

