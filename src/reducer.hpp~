#ifndef UCLUST_HPP_
#define UCLUST_HPP_

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <memory>

// Translate usearch output files into clusters
// clustersIn - input (output of uclust)
// dataFile - output file with cetroids
// clustersOut - output file with clusters
void filter(const std::string& clustersIn, const std::string& dataFile, 
    const std::string& clustersOut);

// Uclust procedure
void procedure(const std::string& centroidsOut, const std::string& clustersOut,
    const char* directory, const int numFiles, float identity = 0.5);

// Divides data in dataFile into parts (files in working directory)
// with defined number of sequences (limit)
// File cannot be larger than 100000000 bytes.
// dataFile - dataset
// directory - working directory
// limit - maximum number of sequences in new files (while dividing)
// iter - internal
// returns number of generated files
int generator(const std::string& dataFile, const char* directory,
    const int limit, const int startId, const int iter);

#endif  // UCLUST_HPP_
