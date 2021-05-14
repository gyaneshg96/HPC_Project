#include<iostream>
#include<string>
#include"utils2.cpp"
#include<vector>
#include<fstream>
#include<cstring>
#include<sstream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include<mpi.h>
#include<cblas.h>

using namespace std;


double *initialize(double *points, long N, int dim, int nclusters, string method){
  
  int mpirank, p;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  double *cluster_centers = (double *)calloc(nclusters*dim, sizeof(double));

  
  if (method == "forgy"){
    #pragma omp parallel for if(nclusters > 1000) 
    for (int i = 0; i < nclusters;i++){
      int val = rand() % N;
      #pragma omp critical
      // memcpy(&cluster_centers[i*dim], &points[val*dim], dim*sizeof(double));
      cblas_dcopy(dim, &points[val*dim], 1, &cluster_centers[i*dim], 1);
      //which clusters
  }
  }
  if (method == "random"){
    int *cluster_nums = (int *)calloc(nclusters, sizeof(int));

    //doubts here
    #pragma omp parallel for schedule(static)
    for (long i = 0; i< N;i++){
      int cluster = rand() % nclusters;
      //add to cluster centers

      //critical section, only allow one thread at a time
      //sllightly erroneous
      #pragma omp critical
      // sum(&cluster_centers[cluster*dim], &points[i*dim], dim);
      cblas_daxpy(dim, 1.0, &points[i*dim], 1, &cluster_centers[cluster*dim], 1);

      cluster_nums[cluster] += 1;
    }

    #pragma omp parallel for if(nclusters > 1000)
    for (int i = 0; i < nclusters; i++)
      // divide(&cluster_centers[i*dim], cluster_nums[i], dim);
      cblas_dscal(dim, 1.0/cluster_nums[i], &cluster_centers[i*dim], 1);
    
    free(cluster_nums);
  }
  if (method == "kmeans++"){

    //first cluster
    int val = rand() % N;
    // cout<<"val "<<val;
    // memcpy(&cluster_centers[0], &points[val*dim], dim*sizeof(double));
    cblas_dcopy(dim, &points[val*dim], 1, &cluster_centers[0], 1);

    // double *maxdist = (double *)calloc(N,sizeof(double));
    for (int i = 1; i < nclusters; i++){
      long farthest = 0;
      double maxdist = 0;
      #pragma omp parallel for 
      for (long j = 0; j < N;j++){
        double mindist = 100000;
        int flag = 0;
        for (int  k = 0; k < i; k++){
          // double currdist = distance(&points[j*dim], &cluster_centers[k*dim], dim);
          // double currdist = distance(&points[j*dim], &cluster_centers[k*dim], dim);
          double currdist = - 2*cblas_ddot(dim, &points[j*dim],1, &cluster_centers[k*dim], 1) +
                            cblas_ddot(dim, &points[j*dim],1, &points[j*dim], 1) +
                            cblas_ddot(dim, &cluster_centers[k*dim],1, &cluster_centers[k*dim], 1);
            currdist = sqrt(currdist);
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
      // memcpy(&cluster_centers[i*dim], &points[farthest*dim], dim*sizeof(double));
      cblas_dcopy(dim, &points[farthest*dim], 1, &cluster_centers[i*dim], 1);
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
  #pragma omp parallel for
  for (int i = 0; i < dim*nclusters;i++){
    for(int j = 0; j < p; j++){
      final_clusters[i] += combine_clusters[i + j*dim*nclusters];  
    }
  }
  // divide(final_clusters, p, dim*nclusters);
  cblas_dscal(dim*nclusters, 1.0/p, final_clusters, 1);
  }

  //send back to other processes
  MPI_Bcast(final_clusters, nclusters*dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(cluster_centers);
  free(combine_clusters);

  return final_clusters;
}

void clusterAssignment(double *points, int* cluster_assignment, double *cluster_centers, int dim, long N, int nclusters){

  #pragma omp parallel for
  for (long point = 0; point < N; point++){

    double min_distance = 100000;
    double min_center = 0;   
    for (int center = 0; center < nclusters; center++){
      // double dist = distance(&cluster_centers[center*dim], &points[point*dim], dim);
         double dist = - 2*cblas_ddot(dim, &points[point*dim],1, &cluster_centers[center*dim], 1) +
                            cblas_ddot(dim, &points[point*dim],1, &points[point*dim], 1) +
                            cblas_ddot(dim, &cluster_centers[center*dim],1, &cluster_centers[center*dim], 1);
          dist = sqrt(dist);
      if (dist < min_distance){
        min_distance = dist;
        min_center = center; 
      }
    }
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

  #pragma omp parallel for
  for (int point = 0; point < N; point++){
    // sum(&temp_centers[cluster_assignment[point]*dim], &points[point*dim], dim);
    cblas_daxpy(dim, 1.0, &points[point*dim], 1, &temp_centers[cluster_assignment[point]*dim], 1);
    cluster_nums[cluster_assignment[point]] += 1;    
  }

  // if (mpirank == 0)
  //   printvec(temp_centers, nclusters*dim);
  
  #pragma omp parallel for if(nclusters > 1000)
  // dont parallelize for small cluster size
  for (int i = 0; i < nclusters; i++){
    // divide(&temp_centers[i*dim], cluster_nums[i], dim);
    cblas_dscal(dim, 1.0/cluster_nums[i], &temp_centers[i*dim], 1);
  }

  double *combine_clusters = NULL;
  if (mpirank == 0){
    combine_clusters = (double *)calloc(nclusters*dim*p, sizeof(double)); 
  }

  //gather all cluster centers to root
  MPI_Gather(temp_centers, nclusters*dim, MPI_DOUBLE, combine_clusters, nclusters*dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (mpirank == 0){
  #pragma omp parallel for
  for (int i = 0; i < dim*nclusters;i++){
    for(int j = 0; j < p; j++){
      cluster_centers[i] += combine_clusters[i + j*dim*nclusters];  
    }
  }
  // divide(cluster_centers, p, dim*nclusters);
  cblas_dscal(dim*nclusters, 1.0/p, cluster_centers, 1);
  }
  MPI_Bcast(cluster_centers, nclusters*dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //shift this to kmeans function
  free(temp_centers);
  free(combine_clusters);
  free(cluster_nums);
}

int hasConverged(int *prev_assignment, int *curr_assignment, long N){
  
  int mpirank, p;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  int final_truth;
  int truth = memcmp(prev_assignment, curr_assignment, N*sizeof(int));
  int* combine_truths = NULL;
  if (mpirank == 0){
    combine_truths = (int *)calloc(p, sizeof(int));
  }

  MPI_Gather(&truth, 1, MPI_INT, combine_truths, 1, MPI_INT, 0, MPI_COMM_WORLD);

  //only stop when all points have converged
  if(mpirank == 0){
    final_truth = 0;
    for(int i = 0; i < p;i++){
      final_truth = final_truth | combine_truths[i]; 
    }
  }

  MPI_Bcast(&final_truth, 1, MPI_INT, 0, MPI_COMM_WORLD);

  free(combine_truths);

  return final_truth;
}

//method which calls KMeans
int trainKMeans(int epochs, double *points, long N, int dim, int nclusters, int *curr_assignment){

  int mpirank, p;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  int epoch = 0;

  //initialize cluster centers
  MPI_Barrier(MPI_COMM_WORLD);

  double tt = MPI_Wtime();
  double *cluster_centers = initialize(points, N, dim, nclusters, "kmeans++");

  if (mpirank == 0){
    // for (int i = 0; i<nclusters; i++)
      // printvec(&cluster_centers[i*dim], dim);
    cout<<"Initialization Time :"<<MPI_Wtime() - tt<<endl;
  }

  clusterAssignment(points, curr_assignment, cluster_centers, dim, N, nclusters);
  
  
  if (mpirank == 0){
    // printvec2(curr_assignment, N);
    cout<<"Epoch\t"<<"Time Taken"<<endl;
  }

  int *prev_assignment = (int *)calloc(N, sizeof(int));
  while(epoch < epochs && hasConverged(prev_assignment, curr_assignment, N)){

    //update the centroids
    // cout<<"BC"<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
    tt = MPI_Wtime();
    updateCentroids(curr_assignment, points, cluster_centers, dim, N, nclusters);
    
    // if (mpirank == 0){
    //   for (int i = 0; i < nclusters; i++)
    //     printvec(&cluster_centers[i*dim], dim);
    // }

    memcpy(prev_assignment, curr_assignment, N*sizeof(int));

    clusterAssignment(points, curr_assignment, cluster_centers, dim, N, nclusters);
    epoch++;

    if (mpirank == 0){
      cout<<epoch<<"\t"<<MPI_Wtime() - tt<<endl;
    }
  }
  // printvec2(curr_assignment,N);
  return epoch;
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
  int nclusters;

  MPI_Barrier(MPI_COMM_WORLD);

  double tt = MPI_Wtime();
  if (argc == 2){
    //use on our dataset
    long n = 20;
    if (n % p == 0){
    for (long i = mpirank; i < 20; i=i+p){
      string filename = "hpcdata/"+ to_string(i*25600)+".csv";
      ifstream f (filename);
      // cout<<mpirank<<" "<<i<<" "<<endl;
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
        if (dim > 0){
            finalarray.insert(finalarray.end(), row.begin(), row.end());
          N++;
        }
    }
    sscanf(argv[1], "%d", &nclusters);
    f.close();
    }
  }
  }
  else {
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
        if (dim > 0){
          if (N % p == mpirank)
            finalarray.insert(finalarray.end(), row.begin(), row.end());
          N++;
        }
    }
    sscanf(argv[2], "%d", &nclusters);
    f.close();
  }

    if (mpirank == 0)
    cout<<"Reading Time: "<<MPI_Wtime() - tt<<endl;
    
    N = finalarray.size();
  
    N = N/dim;
    
    finaldata = finalarray.data();
    // f.close();

    int nthreads = 8;
    omp_set_num_threads(nthreads);
    // if (mpirank == 0)
    // #pragma omp parallel
    // {
    //   cout<<"No of threads "<<omp_get_num_threads()<<" Process "<<mpirank<<endl;
    // }

        
    
    if (mpirank == 0)
      cout<<"Clusters: "<<nclusters<<endl;

    int final_epoch = 1;
    int *curr_assignment = (int *)calloc(N, sizeof(int));
    final_epoch = trainKMeans(100, finaldata, N, dim, nclusters, curr_assignment);
  
    // if (0 == mpirank) {
      // printf("Time elapsed per epoch is %f seconds.\n", elapsed/final_epoch);
    // }
    
  MPI_Finalize();
  return 0;
}