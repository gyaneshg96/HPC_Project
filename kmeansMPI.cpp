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

void output(int* cluster, int size);
void printData(double* data, int size);
void printVec(int* data, int size);
void initializeCentroids(double* dataset, double* centroids, int N, int K, int M);
double getDistance(double* centroid, double* data, int M);
void newCentroids(double* dataset, int* cluster, double* centroids, int* local_cluster_size, int N, int K, int M);
void read_csv(const string& filename, double* dataset, int* label, int N, int M);
double getAccuracy(int* label, int* global_membership, int N);

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);

	if (argc < 5) {
		printf("Usage: ./kmeans <# of clusters> <no of data> <dimension> <filename>\n");
		exit(1);
	}

	int rank, size;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	const int K = atoi(argv[1]);
	const int ROOT = 0;
	const int N = atoi(argv[2]);
	const int M = atoi(argv[3]);
	int elements_per_proc = N / size;

	// Store data as 1D array
	double* dataset = NULL;
	int* label = NULL;

	// 0. read data from data.csv file
	if (rank == ROOT) {
		dataset = (double*) malloc(N * M * sizeof(double));
		label = (int*) malloc(N * sizeof(int));
		read_csv(argv[4], dataset, label, N, M);
		printf("loaded %lu data with dimension of %lu from %s file\n", N, M, "data.csv");
		// printVec(label, N);
	}

	// 1. assign N/P data to each processor
	
	// MPI_Datatype MPI_data;
	// MPI_Type_vector(N, M, M, MPI_DOUBLE, &MPI_data);
	// MPI_Type_commit(&MPI_data);
	
	double* subarray = (double*) malloc(elements_per_proc * M * sizeof(double));
	int* membership = (int*) malloc(elements_per_proc * sizeof(int));
	double* local_means = (double*) malloc(K * M * sizeof(double));
	int* all_memberChanged = (int*) malloc(size * sizeof(int));
	double* all_local_means = (double*) malloc(size * K * M * sizeof(double));
	int* local_cluster_size = (int*) malloc(K * sizeof(int)); // key: cluster id, value: # of data
	int* all_local_cluster_size = (int*) malloc(size * K * sizeof(int));
	
	MPI_Barrier(MPI_COMM_WORLD);
	double tt = MPI_Wtime();

	printf("Rank %d: scatter dataset to each processor\n", rank);
	MPI_Scatter(dataset, elements_per_proc * M, MPI_DOUBLE, subarray, elements_per_proc * M, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	// 2. Node 0 randomly choose K points as cluster means and broadcast
	if (rank == ROOT) initializeCentroids(dataset, local_means, N, K, M);
	
	printf("Rank %d: broadcast k means\n", rank);
	MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Bcast(local_means, K * M, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	int iter = 0;
	double changedPercent = 1.0;


	while(iter < 10000 && changedPercent > 0.0001) {
		// printf("RANK %d - iter%d - local means: ", rank, iter);
		// printData(local_means, K * M);

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
		newCentroids(subarray, membership, local_means, local_cluster_size, elements_per_proc, K, M);

		// 5. globally broadcast all local means for each processor to find the global mean
		MPI_Allgather(local_means, K * M, MPI_DOUBLE, all_local_means, K * M, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(local_cluster_size, K, MPI_INT, all_local_cluster_size, K, MPI_INT, MPI_COMM_WORLD);

		// calculate global means
		// 1.1 calculate sum of means from all processors
		for (int i = 0; i < K * M; ++i) local_means[i] = 0.0; // clear previous local means

		for (int i = 0; i < K; ++i) {
			for (int j = 0; j < M; ++j) {
				for (int p = 0; p < size; ++p) {
					// local mean of cluster i (jth dimension) += local mean of cluster i (jth dimension) on proc p * cluster i sizeon proc p
					local_means[i * M + j] += all_local_means[p * K * M + i * M + j] * all_local_cluster_size[p * K + i];
				}
			}
		}
		
		// 1.2 calculate total size of each cluster
		for (int i = 0; i < K; ++i) local_cluster_size[i] = 0; // clear so can add

		for (int i = 0; i < K; ++i) {
			for (int p = 0; p < size; ++p) {
				local_cluster_size[i] += all_local_cluster_size[p * K + i];
			}
		}

		// 1.3 mean = sum / count
		for (int i = 0; i < K; ++i) {
			for (int j = 0; j < M; ++j) local_means[i * M + j] /= local_cluster_size[i];
		}

		// 6. calculate global member changed
		MPI_Allgather(&memberChanged, 1, MPI_INT, all_memberChanged, 1, MPI_INT, MPI_COMM_WORLD);
		memberChanged = 0;
		for (int p = 0; p < size; ++p) memberChanged += all_memberChanged[p];

		changedPercent = (double) memberChanged / N;

		printf("RANK %d - iteration %d: %d, %7.3f\n", rank, iter, memberChanged, changedPercent);
		// printf("RANK %d - iteration %d:\n", rank, iter);
		// printData(subarray, elements_per_proc);
		// printVec(membership, elements_per_proc);
		// printData(local_means, K);
		iter++;

		MPI_Barrier(MPI_COMM_WORLD); // to make sure iter is in sync. 
	}

	// MPI_Barrier(MPI_COMM_WORLD);

	// Gather Cluster Result
	int* global_membership = NULL;
	if (rank == ROOT) {
		global_membership = (int*) malloc(N * sizeof(int));
	}
	MPI_Gather(membership, elements_per_proc, MPI_INT, global_membership, elements_per_proc, MPI_INT,
				ROOT, MPI_COMM_WORLD);

	if (rank == ROOT) output(global_membership, N);

	tt = MPI_Wtime() - tt;

	if (rank == ROOT) {
		printf("kmeans latency: %e ms\n", tt * 1000);
		printf("kmeans bandwidth: %e GB/s\n", (N * M)/tt/1e9);
		printf("kmeans accuracy: %e%%\n", getAccuracy(label, global_membership, N));
	}
	
	free(dataset);
	free(subarray);
	free(local_means);
	free(all_local_means);
	free(membership);
	free(global_membership);

	MPI_Finalize();

	return 0;
}

void output(int* cluster, int size) {
	for (int i = 0; i < size; ++i) printf("Data %d: cluser %d\n", i, cluster[i]);
}

void printData(double* data, int size) {
	for (int i = 0; i < size; ++i) cout << *data++ << " ";
	cout << endl;
}

void printVec(int* data, int size) {
	for (int i = 0; i < size; ++i) cout << *data++ << " ";
	cout << endl;
}

// random select k centroids
void initializeCentroids(double* dataset, double* centroids, int N, int K, int M) {
	// centroids = (double*) malloc(K * M * sizeof(double));
	srand(time(0));
	int index;
	for (int i = 0; i < K; ++i) {
		index = rand() % N; // randomly choose a data point from dataset
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
void newCentroids(double* dataset, int* cluster, double* centroids, int* local_cluster_size, int N, int K, int M) {
	for (int i = 0; i < K * M; ++i) centroids[i] = 0.0;
	for (int i = 0; i < K; ++i) local_cluster_size[i] = 0;

	int clusterId;

	for (int i = 0; i < N; ++i) {// for each data
		// cluster[i] -> cluster id
		clusterId = cluster[i];
		//printf("\t\tcluster id %d\n", clusterId);
		local_cluster_size[clusterId]++;
		for (int j = 0; j < M; ++j) {// for each dimension
			centroids[clusterId * M + j] += dataset[i * M + j];
		}
	}

	// mean = sum / count
	// assert(newCentroids.size() == k);
	for (int i = 0; i < K; ++i) {
		for (int j = 0; j < M; ++j) {
			if (local_cluster_size[i] > 0) centroids[i*M + j] /= local_cluster_size[i];
		}
	}
}

void read_csv(const string& filename, double* dataset, int* label, int N, int M) {

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
	int row = 0;
    while(getline(myFile, line)) {
        // Create a stringstream of the current line
        stringstream ss(line);
        
        // Extract each integer
        // while(ss >> val){
        //     dataset[i++] = val;
            
        //     // If the next token is a comma, ignore it and move on
        //     if(ss.peek() == ',') ss.ignore();
           
        // }
		while(ss >> val){
			if (i % (M + 1) == 0) label[row] = val;
			else dataset[i - 1 * (row + 1)] = val;
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
           
			++i;
        }

		row++;
    }

    // Close file
    myFile.close();
}

double getAccuracy(int* label, int* global_membership, int N) {
	int correct = 0;
	for (int i = 0; i < N; ++i) correct += label[i] == global_membership[i];
	return (double) correct / N;
}
