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

// typedef vector<double> Data;

void output(const vector<int>& cluster) {
	for (int i = 0; i < cluster.size(); ++i) printf("Data %d: cluser %d\n", i, cluster[i]);
}

void printData(double* data, int size) {
	for (int i = 0; i < size; ++i) cout << *data++ << " ";
	cout << endl;
}

// random select k centroids
void initializeCentroids(double* dataset, double* centroids, int N, int K, int M) {
	centroids = (double*) malloc(K * M * sizeof(double));
	srand(time(0));
	int index;
	for (int i = 0; i < K; ++i) {
		index = rand() % N;
		for (int j = 0; j < M; ++j) centroids[i * M + j] = dataset[index * M + j];
	}
}

double getDistance(double* centroid, double* data, int M) {
	double squareDistance = 0.0;
	for (int i = 0; i < M; ++i) {
		squareDistance += (*data - *centroid) * (*data - *centroid);
		data++;
		centroid++;
	}
	return sqrt(squareDistance);
}


// 1. calculate mean of each cluster as new centroid
// 2. store new centroids
void newCentroids(double* dataset, const vector<int>& cluster, double* centroids, int N, int K, int M) {
	for (int i = 0; i < K * M; ++i) centroids[i] = 0.0;

	vector<int> clusterSize(K, 0);
	int clusterId;

	for (int i = 0; i < N; ++i) {// for each data
		// cluster[i] -> cluster id
		clusterId = cluster[i];
		//printf("\t\tcluster id %d\n", clusterId);
		clusterSize[clusterId]++;
		for (int j = 0; j < M; ++j) {// for each dimension
			centroids[clusterId * M + j] += dataset[i * M + j];
		}
	}

	// mean = sum / count
	// assert(newCentroids.size() == k);
	for (int i = 0; i < K; ++i) {
		for (int j = 0; j < M; ++j)
			centroids[i*M + j] /= clusterSize[i];
	}
}


void read_csv(const string& filename, double* dataset, int N, int M) {

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
	int i = 0;
    while(getline(myFile, line)) {
        // Create a stringstream of the current line
        stringstream ss(line);
        
        // Extract each integer
        while(ss >> val){
            dataset[i++] = val;
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
           
        }
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
	
    const int N = atoi(argv[2]);
    const int M = atoi(argv[3]);
    int elements_per_proc = N / size;

	// Store data as 1D array
	double* dataset;

	// 0. read data from data.csv file
	if (rank == ROOT) {
		dataset = (double*) malloc(N * M * sizeof(double));
		read_csv("data.csv", dataset, N, M);
		printf("loaded %lu data with dimension of %lu from %s file\n", N, M, "data.csv");
	}

	// 1. assign N/P data to each processor
	
	// MPI_Datatype MPI_data;
	// MPI_Type_vector(N, M, M, MPI_DOUBLE, &MPI_data);
	// MPI_Type_commit(&MPI_data);
	
	// vector<Data> subarray(elements_per_proc);
	double* subarray = (double*) malloc(elements_per_proc * M * sizeof(double));
	vector<int> membership(elements_per_proc, -1);
	
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Rank %d: scatter dataset to each processor\n", rank);
	MPI_Scatter(dataset, elements_per_proc * M, MPI_DOUBLE, subarray, elements_per_proc * M, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	// 2. Node 0 randomly choose K points as cluster means and broadcast
	double* local_means = (double*) malloc(K * M * sizeof(double));
	if (rank == ROOT) initializeCentroids(subarray, local_means, elements_per_proc, K, M);
	
	
	printf("Rank %d: broadcast k means\n", rank);
	MPI_Barrier(comm);	
	MPI_Bcast(local_means, K * M, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
	MPI_Barrier(comm);

	int iter = 0;
	double changedPercent = 1.0;

	while(iter < 10000 && changedPercent > 0.001) {
		// 3. for each data point find membership 
		int memberChanged = 0;
		for (int i = 0; i < elements_per_proc; ++i) {
			double minDistance = INT_MAX;
			int closestCentroid;

			for (int j = 0; j < K; ++j) {
				//printf("RANK %d - COMPUTE distance between local means %d and subarray %d\n",rank, j, i);
				//printData(&local_means[j * M], M);
				//printf("RANK %d\n", rank);
				//printData(&subarray[i * M], M);
				double distance = getDistance(&local_means[j * M], &subarray[i * M], M);
				if (distance < minDistance) {
					minDistance = distance;
					closestCentroid = j;
				}
			}

			//printf("Rank %d - Data %d - new centroid: %d\n", rank, i, closestCentroid);

			if (closestCentroid != membership[i]) {
				membership[i] = closestCentroid;
				memberChanged++;
			}
		}

		// 4. Recalculate local means for each cluster in each processor
		newCentroids(subarray, membership, local_means, elements_per_proc, K, M);

		// 5. globally broadcast all local means for each processor to find the global mean
		// vector<Data> all_local_means(size);
		double* all_local_means = (double*) malloc(size * K * M * sizeof(double));
		MPI_Allgather(local_means, K * M, MPI_DOUBLE, all_local_means, K * M, MPI_DOUBLE, MPI_COMM_WORLD);

		// calculate global means
		// 1.1 calculate sum of means from all processors
		// vector<Data> global_means(K, vector<double>(M, 0.0));
		for (int i = 0; i < K * M; ++i) local_means[i] = 0.0; // clear previous local means

		
		for (int i = 0; i < K; ++i) {
			for (int j = 0; j < M; ++j) {
				for (int p = 0; p < size; ++p) {
					local_means[i * M + j] += all_local_means[p * K * M + i * M + j];
				}
			}
		}
		

		// 1.2 mean = sum of means / number of processors
		for (int i = 0; i < K; ++i) {
			for (int j = 0; j < M; ++j) local_means[i * M + j] /= size;
		}

		changedPercent = (double) memberChanged / elements_per_proc;
		printf("iteration %d: %d, %7.3f\n", iter, memberChanged, changedPercent);
		iter++;
	}

	// MPI_Type_free(&MPI_data);
	MPI_Finalize();

	return 0;
}
