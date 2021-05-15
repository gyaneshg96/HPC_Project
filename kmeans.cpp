#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <climits>

using namespace std;

typedef vector<double> Data;

/****** K-means algorithm
random select k data centroids
while there is change in k centroids:

	for each data in dataset:
		mini distacnce
		for each centroid in k:
			cal euclidean distance
			update minimum distance

		assign data to nearest cluster

	calculate the mean of each cluster as next centroid (each cluster stays on one thread)
*******/

class KMeans {
public:
	KMeans(int k_in) : k(k_in) {}

	void solution(const vector<Data>& dataset) {
		cluster.resize(dataset.size(), -1);
		initializeCentroids(dataset);
		// printf("created %lu centroids\n", centroids.size());
		// for (int i = 0; i < centroids.size(); ++i) {
		// 	printf("Centroid %d:", i);
		// 	printData(centroids[i]);
		// }

		int iter = 0; // for debug use
		double changedPercent = 1.0;

		while(iter < 10000 && changedPercent >= 0.001) {
			int memberChanged = 0; // count number of data that has changed memebership

			for (int i = 0; i < dataset.size(); ++i) {
				double minDistance = INT_MAX;
				int closestCentroid;

				for (int j = 0; j < centroids.size(); ++j) {
					double distance = getDistance(centroids[j], dataset[i]);
					if (distance <= minDistance) {
						minDistance = distance;
						closestCentroid = j;
					}
				}

				// printf("new centroid: %d\n", closestCentroid);

				if (closestCentroid != cluster[i]) {
					cluster[i] = closestCentroid;
					memberChanged++;
				}
				
				// printf("Data %d belongs to cluster %d\n", i, closestCentroid);
			}

			// calculate the mean of each cluster as next centroid (each cluster stays on one thread)
			newCentroids(dataset);
			
			changedPercent = (double) memberChanged / dataset.size();
			printf("iteration %d: %d, %7.3f\n", iter, memberChanged, changedPercent);
			iter++;
		}
	}

	void output() {
		for (int i = 0; i < cluster.size(); ++i) printf("Data %d: cluser %d\n", i, cluster[i]);
	}
private:

	void printData(const Data& data) {
		for (int i = 0; i < data.size(); ++i) cout << data[i] << " ";
		cout << endl;
	}

	// random select k centroids
	void initializeCentroids(const vector<Data>& dataset) {
		srand(time(0));
		int index;
		for (int i = 0; i < k; ++i) {
			index = rand() % dataset.size();
			centroids.push_back(dataset[index]);
		}
	}

	double getDistance(const Data& centroid, const Data& data) {
		assert(centroid.size() == data.size());
		double squareDistance = 0.0;
		for (int i = 0; i < data.size(); ++i) {
			squareDistance += (data[i] - centroid[i]) * (data[i] - centroid[i]);
		}
		return sqrt(squareDistance);
	}


	// 1. calculate mean of each cluster as new centroid
	// 2. store new centroids
	void newCentroids(const vector<Data>& dataset) {
		vector<Data> newCentroids(centroids.size(), vector<double>(dataset[0].size(), 0.0));
		vector<int> clusterSize(centroids.size(), 0);
		int clusterId;

		for (int i = 0; i < dataset.size(); ++i) {// for each data
			// cluster[i] -> cluster id
			clusterId = cluster[i];
			clusterSize[clusterId]++;
			for (int j = 0; j < dataset[i].size(); ++j) {// for each dimension
				newCentroids[clusterId][j] += dataset[i][j];
			}
		}

		// mean = sum / count
		assert(newCentroids.size() == k);
		for (int i = 0; i < k; ++i) {
			for (int j = 0; j < newCentroids[i].size(); ++j)
				newCentroids[i][j] /= clusterSize[i];
		}

		centroids = newCentroids;
	}

	int k;
	vector<Data> centroids;
	vector<int> cluster;
};


void read_csv(const string& filename, vector<Data>& dataset) {

    // Create an input filestream
    ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) {
    	printf("Could not open file\n");
    	exit(1);
    }

    // Helper vars
    string line;
    int val;

    // Read data, line by line
    while(getline(myFile, line)) {
        // Create a stringstream of the current line
        stringstream ss(line);
        
        Data data;
        
        // Extract each integer
        while(ss >> val){
            
            data.push_back(val);
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
           
        }
        dataset.push_back(data);
    }

    // Close file
    myFile.close();
}


int main(int argc, char* argv[]) {
	// read data from data.csv file
	vector<Data> dataset;
	read_csv("data.csv", dataset);
	printf("loaded %lu data with dimension of %lu from %s file\n", dataset.size(), dataset[0].size(), "data.csv");

	// test
	int k = 3;
	KMeans kmeans(k);

	kmeans.solution(dataset);
	kmeans.output();

	return 0;
}