#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <time.h>
#include <ctime>
#include <getopt.h>

#include "reducer.hpp"

int main(int argc, char **argv) {
  // FLAGS
  int c;
  std::string input;
  std::string centFile = "centroids.txt";
  std::string clustersFile = "clusters.txt";
  float identity = 0.5;
  bool flagIn = 0;
  bool flagCent = 0;
  bool flagClu = 0;
    
  while (1) {
    int this_option_optind = optind ? optind : 1;
        
    int option_index = 0;
    static struct option long_options[] = {
        {"in",  required_argument, 0, 'i'},
        {"dataset",     required_argument, 0,  'd' },
        {"clusters",  required_argument, 0,  'c' },   
        {"identity",    required_argument, 0,  'p' },
        {0,         0,                 0,  0 }
    };

    c = getopt_long(argc, argv, "i:o:c:p:",
                 long_options, &option_index);
        
    if (c == -1)
        break;

    switch (c) {
      case 'i':
        input.assign(optarg);
        flagIn = 1;
        break;
      case 'd':
       	centFile.assign(optarg);
       	flagCent = 1;
       	break;
      case 'c':
        clustersFile.assign(optarg);
        flagClu = 1;
        break;
      case 'p':
        identity = atof(optarg);
        break;
      case '?':
        break;
      }
  }

  if ((flagIn & flagCent & flagCent) == 0) {
    std::cout << "Please provide input sequence (--in), dataset" 
        << " (--dataset) and clusters (--clusters) ";
	exit(0);
  }

  // Create working directory
  char temp[] = "tmpXXXXXX";
  char *f = mkdtemp(temp);
  if (f == NULL) {
  	std::cout << "Failed";
  	exit(0);
  }

  // Searching for a sequence
  procedure(input, centFile, clustersFile, identity, temp);

  // Deleting working directory and files
  std::ostringstream oss;
  oss << "rm " << temp << "/*";
  system(oss.str().c_str());

  oss.str("");
  oss.clear();
  oss << "rmdir " << temp;
  system(oss.str().c_str());
  return 0;
}

void procedure(const std::string& input, const std::string& centFile, 
	const std::string& clustersFile, float identity, const char* temp) {

  // Read sequence from input file, one by one
  std::ifstream fileData(input);
  std::string line;
  std::ostringstream dataSeq;
  std::string seqName;

  std::ostringstream oss;
  oss << temp << "/datasetFile.txt";
  std::string dataSetFile = oss.str();
	
  getline(fileData, line);
  dataSeq << line << '\n';
  seqName = line.substr(1, line.size());

  int number;
  bool runFlag = 0;
  bool endFlag = 0;

  int counter = 0;
  double totalTime = 0;
  while(true) {
    if (getline(fileData, line) == 0) {
		runFlag = 1;
		endFlag = 1;
	} else if (line.compare(0, 1, ">") == 0) {
		runFlag = 1;
	}

	if (runFlag) {
        std::cout << seqName << ":\n";
		oss.str("");
		oss << temp << "/seq.fasta";
		std::ofstream myfile(oss.str());
		
		myfile << dataSeq.str();
		myfile.close();

		dataSeq.str("");
		dataSeq.clear();

		clock_t begin = clock();

		// Find potential clusters
		std::vector<std::string> potClusters;

		selectSeq(oss.str(), centFile, potClusters, temp, identity);

		for (std::vector<std::string>::iterator it = potClusters.begin(); 
									it != potClusters.end(); it++) {
	        std::vector<std::string> results;
			number = getDataset(clustersFile, *it, dataSetFile);

			if (number < 2) {
			    std::cout << "  " << *it << "\n";
				results.push_back(*it);

				if (seqName.compare(*it) == 0) {
					counter++;
				}
				continue;
			}
							
			selectSeq(oss.str(), dataSetFile, results, temp, identity);

			for (std::vector<std::string>::iterator it_in = results.begin(); 
									it_in != results.end(); it_in++) {
				std::cout << "  " << *it_in << "\n";
				if (seqName.compare(*it_in) == 0)
					counter++;
			}
		} 

		totalTime += (1.0 * (clock() - begin)) / CLOCKS_PER_SEC;

		if (endFlag) {
			break;
		} else {
			seqName = line.substr(1, line.size());
		}
	}

	runFlag = 0;

	dataSeq << line << '\n';
  }

  std::cout << "Total time: " << totalTime << "\n";
  fileData.close();
  std::cout <<"Correct: " <<counter; 

  return;
}

int getDataset(const std::string& fileIn, const std::string& clusterName, 
  const std::string& out) {
  std::ifstream file(fileIn);
  std::ofstream fileOut(out);

  int number = 0;

  std::string lineComm = "Cluster:\t" + clusterName;
  
  std::string line;
  while(getline(file, line)) {
	if ((line.compare(0, lineComm.size(), lineComm) == 0) && 
	  (line.size() == lineComm.size())) {

	  getline(file, line);
	  while (!line.empty()) {
		if (line.compare(0, 1, ">") == 0) number++;
        fileOut << line << "\n";
		getline(file, line);
	  }

	  break;
	} else {
	  while (!line.empty()) {
		getline(file, line);
	  }
	}
  }

  fileOut.close();
  file.close();
  return number;
}

void selectSeq(const std::string& fileIn, const std::string& fileCentroids, 
    std::vector<std::string>& potClusters, const char* directory, float identity) {
  // Global alignment to cluster
  std::ostringstream command;
  command << "./usearch -usearch_global " << fileIn << " -db " << fileCentroids
				<< " -gapopen 5.0I/1.0E -gapext 2.0I/0.5E -id "<< identity 
				<<" -maxaccepts 1  -blast6out " <<directory <<"/results.txt -strand plus --quiet";

  system(command.str().c_str());

  // take top 5 
  std::string line;

  std::ostringstream oss;
  oss << directory << "/results.txt";
  std::ifstream file(oss.str());
  std::string value;

  for (int i = 0; i < 3 && getline(file, line); i++) {
	std::istringstream iss;
    iss.str(line);
    iss >> value >> value;
    potClusters.push_back(value);
  }
  file.close();

  command.str("");
  command << "rm " << directory << "/results.txt";
  system(command.str().c_str());
}
