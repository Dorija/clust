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
#include <getopt.h>

#include "reducer.hpp"

int main(int argc, char** argv) {
  //OPTIONS
  int c;
  std::string input;
  std::string cent = "centroids.txt";
  std::string clusters = "clusters.txt";
  float identity = 0.5;
  bool flag = 0;

  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
        {"in",  required_argument, 0, 'i'},
        {"out",     required_argument, 0,  'o' },
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
                flag = 1;
            break;
       	case 'o':
       	    cent.assign(optarg);
       	    break;
        case 'c':
            clusters.assign(optarg);
            break;
        case 'p':
            identity = atof(optarg);
            break;
        case '?':
            break;
    }
  }

  if (flag == 0) {
    std::cout << "Please provide input dataset (--in) ";
	exit(0);
  }

  // creating temporary folder
  char temp[] = "tmpXXXXXX";
  char *f = mkdtemp(temp);
  
  if (f == NULL) {
	std::cout << "Failed";
  	exit(0);
  }

  // Generate data in temp directory
  int numFiles = generator(input, temp, 50000, 0, 0);

  // Create clusters from data in temp directory
  procedure(cent, clusters, temp, numFiles, identity);

  //delete temp directory;
  std::ostringstream command;
  command << "rm " << temp << "/*";
  system(command.str().c_str());
  command.str("");

  command << "rmdir " << temp;
  system(command.str().c_str());
  command.clear();

  return 0;
}

int generator(const std::string& dataFile, const char *directory, const int limit, const int startId, const int iter) {
  using namespace std;

  ifstream data(dataFile);
  string line;
  int counter = 0;
  int fileId = startId;

  bool flag = 0;

  ostringstream dataStream;

  long maxSize = 10000000L;
  long sizeData;
  while(getline(data,line)) {
    if (line.compare(0, 1, ">") == 0) {
    	counter++;

		if (counter == limit) {
			sizeData = dataStream.tellp();

			//if data is to big, divide into partitions
			if (sizeData > maxSize) {
				flag = 1;
				ostringstream spec;
				spec << directory << "/tmp" << iter;
				ofstream file(spec.str(), std::ofstream::app);
					
				file << dataStream.str();

				file.close();
				dataStream.clear();
				dataStream.str("");
				counter = 0;

			} else {
				ostringstream spec;
				spec << directory << "/data" << fileId;
				ofstream file(spec.str());

				file << dataStream.str();

				file.close();
				dataStream.clear();
				dataStream.str("");
				counter = 0;
				fileId++;
			}	
		}
	}

    dataStream << line << "\n";
  }

  data.close();

  if (counter > 0) {
  	sizeData = dataStream.tellp();

	//if data is to big, divide into partitions
	if (sizeData > maxSize) {
		flag = 1;
		ostringstream spec;
		spec << directory << "/tmp" << iter;
		ofstream file(spec.str(), std::ofstream::app);
					
		file << dataStream.str();

	    file.close();
	    dataStream.clear();
	    dataStream.str("");

	} else {
		ostringstream spec;
		spec << directory << "/data" << fileId;
		ofstream file(spec.str());

		file << dataStream.str();

		file.close();
        dataStream.clear();
		dataStream.str("");
		fileId++;
	}				
  }

  if (flag == 1) {
	ostringstream fileName;
	fileName << directory << "/tmp" << iter;
	if (iter < 3) {
		fileId = generator(fileName.str(), directory, limit/2, fileId, iter + 1);
	} else { 
		fileId = generator(fileName.str(), directory, limit - 0.1 * limit, fileId, iter + 1);
	}
  }
  return fileId;
}

void procedure(const std::string& centroidsOut, const std::string& clustersOut, const char *directory, const int numFiles, const float identity) {
  using namespace std;

  //paralel proccesing
  #pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < numFiles; i++) {
	ostringstream dataF;
	ostringstream centroidF;
	ostringstream clusterF;
	ostringstream outF;

	dataF << directory <<"/data" << i;
	centroidF << directory <<"/cent" << i;
	clusterF << directory <<"/clusters" << i;
	outF << directory <<"/out" << i;
		
	ostringstream command;
	
	command << "./usearch -sort length -gapopen 5.0I/1.0E -gapext 2.0I/0.5E -id " 
		<< identity << " --maxrejects 8 -cluster_fast " << dataF.str() 
		<< " -centroids " << centroidF.str() << " -uc "<< clusterF.str() << " --quiet";
	int status = system(command.str().c_str());

	filter(clusterF.str(), dataF.str(),  outF.str());
  }

  //clear output files
  ofstream fileT(clustersOut);
  fileT.close();
  ofstream fileL(centroidsOut);
  fileL.close();

  //conect all file in one out file
  ostringstream oss;
  oss << "cat " << directory << "/out* >> " << clustersOut;
  system (oss.str().c_str());
  oss.str("");
  oss << "cat " << directory << "/cent* >> " << centroidsOut;
  system (oss.str().c_str());
}

void filter(const std::string& clustersIn, const std::string& dataFile, const std::string& clustersOut) {
  using namespace std;
  //reading UCLUST-formats to identify clusters
  ifstream fileIn(clustersIn);
  string line;
  istringstream iss;
  string value;
  string centroid;
  map <string, int> dict;
  map<int, int> clusters;
  vector<int> centroids;

  int identifier = 1;
  while (getline(fileIn, line)) {
    if (line.compare(0, 1, "H") == 0) {
		iss.str(line);

		//skip first 8 columms
		for (int i = 0; i < 8; i++) iss >> value;
		  iss >> value; //sequence
		  iss >> centroid; //cluster

		  if (dict[value] == 0) 
			dict[value] = identifier++;
				
		  if (dict[centroid] == 0)
			dict[centroid] = identifier++;

          clusters.insert(pair<int, int>(dict[value], dict[centroid]));

		  iss.clear();

	} else if (line.compare(0, 1, "S") == 0) {
		iss.str(line);

		//skip first 8 columms
		for (int i = 0; i < 8; i++) iss >> value;

		iss >> centroid; //cluster

		if (dict[centroid] == 0)
			dict[centroid] = identifier++;

		clusters.insert(pair<int, int>(dict[centroid], dict[centroid]));
		
		iss.clear();
		centroids.push_back(dict[centroid]);
	}
  }

  fileIn.close();
  map<int, std::shared_ptr<std::ostringstream>> stream_map;

  for (vector<int>::iterator it = centroids.begin(); it != centroids.end(); it++) {
		stream_map.insert(pair<int, shared_ptr<ostringstream>>(*it, shared_ptr<ostringstream>(new std::ostringstream())));
  }

  //reading dataset to divide them into clusters
  ifstream data(dataFile);
  int centroidID = 0;
  while(getline(data, line)) {
	if (line.compare(0, 1, ">") == 0) {
		value = line.substr(1, line.size());

		//which cluster it belongs
		centroidID = clusters[dict[value]];
	}

	(*stream_map[centroidID]) << line << "\n";
  }

  data.close();

  //make reversed dict
  map<int, const string*> dictReverse;
  for (map<string, int>::iterator it = dict.begin(); it != dict.end(); it++) {
	dictReverse.insert(pair<int, const string*>(it->second, &(it->first)));
  }

  //write into outputfile
  ofstream fileOut(clustersOut, std::fstream::app);
  for (map<int, shared_ptr<ostringstream>>::iterator it = stream_map.begin(); it != stream_map.end(); it++) {
	fileOut << "Cluster:\t" << *(dictReverse[it->first]) << "\n";
	fileOut << (*(it->second)).str();
	fileOut << "\n";
  }

  fileOut.close();
  stream_map.clear();
}
