#include<iostream>
#include<string>
#include"utils.cpp"
#include<vector>
#include<fstream>
#include<cstring>
#include<sstream>
#include<mpi.h>

using namespace std;



double *initialize(double *points, long N, int dim, int nclusters, string method){
  // double **cluster_centers = (double **)malloc(nclusters*sizeof(double *));


  int mpirank, p;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  double *cluster_centers = (double *)calloc(nclusters*dim, sizeof(double));

  
  /*
  for (int i = 0; i < nclusters; i++){
    cluster_centers[i] = (double *)calloc(dim, sizeof(double));
  }*/

  if (method == "forgy"){
    for (int i = 0; i < nclusters;i++){
      int val = rand() % N;
      // memcpy(cluster_centers[i], points[val], dim*sizeof(double));
      memcpy(&cluster_centers[i*dim], &points[val*dim], dim*sizeof(double));
      //which clusters
  }
  }
  if (method == "random"){
    int *cluster_nums = (int *)calloc(nclusters, sizeof(int));
    for (long i = 0; i< N;i++){
      int cluster = rand() % nclusters;
      //add to cluster centers
      sum(&cluster_centers[cluster*dim], &points[i*dim], dim);
      cluster_nums[cluster] += 1;
    }
    for (int i = 0; i < nclusters; i++)
      // cluster_centers[i] /= cluster_nums[i];
      divide(&cluster_centers[i*dim], cluster_nums[i], dim);
    
    free(cluster_nums);
  }
  if (method == "kmeans++"){

    //first cluster
    int val = rand() % N;
    // cout<<"val "<<val;
    memcpy(&cluster_centers[0], &points[val*dim], dim*sizeof(double));

    // double *maxdist = (double *)calloc(N,sizeof(double));
    for (int i = 1; i < nclusters; i++){
      long farthest = 0;
      double maxdist = 0;
      for (long j = 0; j < N;j++){
        double mindist = 100000;
        int flag = 0;
        for (int  k = 0; k < i; k++){
          double currdist = distance(&points[j*dim], &cluster_centers[k*dim], dim);
          if (currdist > 0.00001){ //do not choose the point again
            mindist = min(mindist, currdist);
          }
          else{
            flag = 1;
            //not a candidate
          }
        }
        if (flag == 1)
          continue;
        //for now we are just taking maximum
        if (maxdist < mindist){
          maxdist = mindist;
          farthest = j;
          // cout<<mindist<<" "<<farthest<<endl;
        }
      }
      // cout<<farthest<<endl;
      memcpy(&cluster_centers[i*dim], &points[farthest*dim], dim*sizeof(double));
    }
  }

  //each process has a copy of cluster centers
  //take mean of each and distribute back to the processes

  
  double *final_clusters = (double*)calloc(nclusters*dim, sizeof(double));
  double *combine_clusters = NULL;
  if (mpirank == 0){
    combine_clusters = (double *)calloc(nclusters*dim*p, sizeof(double)); 
  }

  //gather all cluster centers to root
  
  MPI_Gather(cluster_centers, nclusters*dim, MPI_DOUBLE, combine_clusters, nclusters*dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if (mpirank == 0){
  for (int i = 0; i < dim*nclusters;i++){
    for(int j = 0; j < p; j++){
      final_clusters[i] += combine_clusters[i + j*dim*nclusters];  
    }
  }
  divide(final_clusters, p, dim*nclusters);
  }

  //send back to other processes
  MPI_Bcast(final_clusters, nclusters*dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(cluster_centers);
  free(combine_clusters);

  return final_clusters;
}

void clusterAssignment(double *points, int* cluster_assignment, double *cluster_centers, int dim, long N, int nclusters){

  for (long point = 0; point < N; point++){

    double min_distance = 100000;
    double min_center = 0;   

    for (int center = 0; center < nclusters; center++){
      double dist = distance(&cluster_centers[center*dim], &points[point*dim], dim);
      if (dist < min_distance){
        min_distance = dist;
        min_center = center; 
      }
    }
    // cout<<min_center<<endl;
    cluster_assignment[point] = min_center;
  }
}
void updateCentroids(int *cluster_assignment, double *points, double *cluster_centers, int dim, long N, int nclusters){
  
  int mpirank, p;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  memset(cluster_centers, 0, nclusters*dim*sizeof(double));
  double* temp_centers = (double*)calloc(nclusters*dim, sizeof(double));

  int *cluster_nums = (int *)calloc(nclusters, sizeof(int));
  for (int point = 0; point < N; point++){
    sum(&temp_centers[cluster_assignment[point]*dim], &points[point*dim], dim);
    cluster_nums[cluster_assignment[point]] += 1;    
  }
  
  for (int i = 0; i < nclusters; i++){
    divide(&temp_centers[i*dim], cluster_nums[i], dim);
  }

  double *combine_clusters = NULL;
  if (mpirank == 0){
    combine_clusters = (double *)calloc(nclusters*dim*p, sizeof(double)); 
  }

  //gather all cluster centers to root
  MPI_Gather(temp_centers, nclusters*dim, MPI_DOUBLE, combine_clusters, nclusters*dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (mpirank == 0){
  for (int i = 0; i < dim*nclusters;i++){
    for(int j = 0; j < p; j++){
      cluster_centers[i] += combine_clusters[i + j*dim*nclusters];  
    }
  }
  divide(cluster_centers, p, dim*nclusters);
  }
  MPI_Bcast(cluster_centers, nclusters*dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // cout<<mpirank<<" !! "<<endl;
  //shift this to kmeans function
  free(temp_centers);
  free(combine_clusters);
  free(cluster_nums);
}

//method which calls KMeans
void trainKMeans(int epochs, double *points, long N, int dim, int nclusters){

  int mpirank, p;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  int epoch = 0;

  //initialize cluster centers
  double *cluster_centers = initialize(points, N, dim, nclusters, "kmeans++");

  if (mpirank == 0)
    cout<<"Init"<<endl;

  int *curr_assignment = (int *)calloc(N, sizeof(int));


  clusterAssignment(points, curr_assignment, cluster_centers, dim, N, nclusters);
  
  if (mpirank == 0){
  for (int i = 0; i < nclusters; i++)
    printvec(cluster_centers + i*dim, dim);
  
  cout<<"Cluster Assigned"<<endl;
  }

  int *prev_assignment = (int *)calloc(N, sizeof(int));
  while(epoch < epochs && memcmp(prev_assignment, curr_assignment, N*sizeof(int)) ){

    //update the centroids
    updateCentroids(curr_assignment, points, cluster_centers, dim, N, nclusters);
    
    if (mpirank == 0){
      for (int i = 0; i < nclusters; i++)
        printvec(&cluster_centers[i*dim], dim);
    }

    memcpy(prev_assignment, curr_assignment, N*sizeof(int));

    clusterAssignment(points, curr_assignment, cluster_centers, dim, N, nclusters);
    epoch++;

    if (mpirank == 0){
    // printvec2(curr_assignment, N);
    cout<<"Epoch "<<epoch<<endl;
    }
  }
  // printvec2(curr_assignment,N);
  // return curr_assignment;
}

int main(int argc, char* argv[]){
  double* finaldata;
  string line;
  vector<double> finalarray;

  int mpirank, p;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  int dim = 0;
  long N = 0;

  if (argc == 1){
    //random generate
  }
  else {
    if (mpirank == 0){
      ifstream f (argv[1]);   /* open file */
    if (!f.is_open()) {     /* validate file open for reading */
        perror (("error while opening file " + string(argv[1])).c_str());
        return 1;
    }

    while (getline(f, line)) {         /* read each line */
        string val;                     /* string to hold value */
        stringstream s (line);
        vector<double> row;
        while (getline(s, val, ',')){
            try {
              double dval = stod(val);
              row.push_back(dval);
            }
            catch (exception& e){
              //do nothing
            }
          }
        if (dim == 0)
        dim = row.size();
        if (dim > 0)
          finalarray.insert(finalarray.end(), row.begin(), row.end());
    }

    N = finalarray.size();

    N = N/dim;
    
    finaldata = finalarray.data();
    f.close();
    }

    MPI_Bcast(&N, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double *eachdata = (double *)calloc(dim*N/p, sizeof(double));

    MPI_Scatter(finaldata, N*dim/p, MPI_DOUBLE, eachdata, N*dim/p, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    /*
    for (int i = 0; i < N*dim/p; i++){
      cout<<eachdata[i]<<" "<<i<<endl;
    }*/
    
    // int nclusters = sscanf(argv[2]);
    
    int nclusters = 3;

    //we are sending N/p once and never again
    // int* assignment = trainKMeans(100, eachdata, N/p, dim, nclusters);
    trainKMeans(100, eachdata, N/p, dim, nclusters);
  }
  MPI_Finalize();
  return 0;
}