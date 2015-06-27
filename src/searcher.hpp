#ifndef DATA_HPP_
#define DATA_HPP_

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>


// The function chose top 3 clusters represented by the centroids
// which match the best with the sequence
// fileIn - file contrains the sequence
// fileCentroids - file constains all centroids from the database
// potClusters - list of potential clusters
// directory - working directory
// identity - identity value
void selectSeq(const std::string& fileIn, const std::string& fileCentroids,
    std::vector<std::string>& potClusters, const char* directory, float identity);

// The function creates a file in FASTA format, which consists of the sequences 
// from the wanted cluster 
// fileIn - file constais clusters 
// clusterName - the name of the wanted cluster
// out - the name of the output file
// returns number of sequences in the cluster
int getDataset(const std::string& fileIn, const std::string& clusterName,
    const std::string& out);

// The procedure runs an identification of the sequences given in the file input
// input - fasta format with input sequences
// centFile - dataset file in FASTA format
// clustersFile - file which contrains a clusters of the dataset
// identity - identity for alignment
// temp - internal
void procedure(const std::string& input, const std::string& centFile,
	const std::string& clustersFile, float identity, const char* temp);

#endif  // DATA_HPP_