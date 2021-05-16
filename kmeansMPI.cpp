#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <climits>
#include <assert.h>
#include <mpi.h>

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


/*void solution(const vector<Data>& dataset) {
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
}*/

void output(const vector<Data>& centroids, const vector<int>& cluster) {
	for (int i = 0; i < cluster.size(); ++i) printf("Data %d: cluser %d\n", i, cluster[i]);
}

void printData(const Data& data) {
	for (int i = 0; i < data.size(); ++i) cout << data[i] << " ";
	cout << endl;
}

// random select k centroids
void initializeCentroids(const vector<Data>& dataset, vector<Data>& centroids) {
	srand(time(0));
	int index;
	for (int i = 0; i < centroids.size(); ++i) {
		index = rand() % dataset.size();
		centroids[i] = dataset[index];
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
void newCentroids(const vector<Data>& dataset, vector<int>& cluster, vector<Data>& centroids) {
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
	// assert(newCentroids.size() == k);
	for (int i = 0; i < centroids.size(); ++i) {
		for (int j = 0; j < newCentroids[i].size(); ++j)
			newCentroids[i][j] /= clusterSize[i];
	}

	centroids = newCentroids;
}


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
	MPI_Init(&argc, &argv);

	if (argc < 4) {
		printf("Usage: ./kmeans <# of clusters> <no of data> <dimension>\n");
		exit(1);
	}

	const int K = atoi(argv[1]);
	const int ROOT = 0;

	int rank, size;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	
	vector<Data> dataset;
        const int N = atoi(argv[2]);
	const int M = atoi(argv[3]);
        int elements_per_proc = N / size;

	// 0. read data from data.csv file
	if (rank == ROOT) {
		read_csv("data.csv", dataset);
		printf("loaded %lu data with dimension of %lu from %s file\n", dataset.size(), dataset[0].size(), "data.csv");
	}

	// 1. assign N/P data to each processor
	//const int N = dataset.size();
	//const int M = dataset[0].size();
	//int elements_per_proc = dataset.size() / size;
	
	MPI_Datatype MPI_data;
	MPI_Type_vector(N, M, M, MPI_DOUBLE, &MPI_data);
	MPI_Type_commit(&MPI_data);
	
	vector<Data> subarray(elements_per_proc);
	vector<int> membership(elements_per_proc, -1);
	
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Rank %d: scatter dataset to each processor\n", rank);
	MPI_Scatter(&dataset[0], elements_per_proc, MPI_data, &subarray[0], elements_per_proc, MPI_data, ROOT, MPI_COMM_WORLD);
	/*for (int i = 0; i < elements_per_proc; ++i) {
		subarray[i] = dataset[rank * elements_per_proc + i];
	}*/

	// 2. Node 0 randomly choose K points as cluster means and broadcast
	vector<Data> local_means(K);
	if (rank == ROOT) initializeCentroids(subarray, local_means);
	
	
	printf("Rank %d: broadcast k means\n", rank);
	MPI_Barrier(comm);	
	MPI_Bcast(&local_means[0], K, MPI_data, ROOT, MPI_COMM_WORLD);
	MPI_Barrier(comm);

	int iter = 0;
	double changedPercent = 1.0;

	while(iter < 10000 && changedPercent > 0.001) {
		// 3. for each data point find membership 
		int memberChanged = 0;
		for (int i = 0; i < subarray.size(); ++i) {
			double minDistance = INT_MAX;
			int closestCentroid;

			for (int j = 0; j < local_means.size(); ++j) {
				printf("RANK %d - COMPUTE distance between local means %d and subarray %d\n",rank, j, i);
				printData(local_means[j]);
				printf("RANK %d\n", rank);
				printData(subarray[i]);
				double distance = getDistance(local_means[j], subarray[i]);
				if (distance <= minDistance) {
					minDistance = distance;
					closestCentroid = j;
				}
			}

			// printf("new centroid: %d\n", closestCentroid);

			if (closestCentroid != membership[i]) {
				membership[i] = closestCentroid;
				memberChanged++;
			}
		}

		// 4. Recalculate local means for each cluster in each processor
		newCentroids(subarray, membership, local_means);

		// 5. globally broadcast all local means for each processor to find the global mean
		vector<Data> all_local_means(size * K);
		MPI_Allgather(&local_means[0], K, MPI_data, &all_local_means[0], K, MPI_data, MPI_COMM_WORLD);

		// calculate global means
		// 1.1 calculate sum of means from all processors
		vector<Data> global_means(K, vector<double>(M, 0.0));
		for (int i = 0; i < all_local_means.size(); ++i) {
			for (int j = 0; j < M; ++j) {
				global_means[i % K][j] += all_local_means[i][j];
			}
		}

		// 1.2 mean = sum of means / number of processors
		for (int i = 0; i < K; ++i) {
			for (int j = 0; j < M; ++j) global_means[i][j] /= size;
		}

		local_means = global_means;
		
	}

	MPI_Type_free(&MPI_data);
	MPI_Finalize();

	return 0;
}
