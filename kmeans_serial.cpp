#include<iostream>
#include<string>
#include"utils2.cpp"
#include<vector>
#include<fstream>
#include<cstring>
#include<sstream>
#include "utils.h"

using namespace std;



double *initialize(double *points, long N, int dim, int nclusters, string method){
  // double **cluster_centers = (double **)malloc(nclusters*sizeof(double *));
  double *cluster_centers = (double *)malloc(nclusters*dim*sizeof(double));
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
  return cluster_centers;
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
    cluster_assignment[point] = min_center;
  }
}
void updateCentroids(int *cluster_assignment, double *points, double *cluster_centers, int dim, long N, int nclusters){
  
  memset(cluster_centers, 0, nclusters*dim*sizeof(double));

  int *cluster_nums = (int *)calloc(nclusters, sizeof(int));
  for (int point = 0; point < N; point++){
    sum(&cluster_centers[cluster_assignment[point]*dim], &points[point*dim], dim);
    cluster_nums[cluster_assignment[point]] += 1;    
  }
  // printvec(cluster_centers[0], dim);

  for (int i = 0; i < nclusters; i++){
    divide(&cluster_centers[i*dim], cluster_nums[i], dim*nclusters);
  }

  free(cluster_nums);
}


int trainKMeans(int epochs, double *points, long N, int dim, int nclusters, int* curr_assignment){  
  int epoch = 0;
  double *cluster_centers = initialize(points, N, dim, nclusters, "kmeans++");

  cout<<"Init"<<endl;

  clusterAssignment(points, curr_assignment, cluster_centers, dim, N, nclusters);
  
  // for (int i = 0; i < nclusters; i++)
  //   printvec(cluster_centers + i*dim, dim);

  cout<<"Cluster Assigned"<<endl;
  //all zero
  int *prev_assignment = (int *)calloc(N, sizeof(int));

  // while(compareAssignment(prev_assignment, curr_assignment, N) && epoch < epochs){
  while(epoch < epochs && memcmp(prev_assignment, curr_assignment, N*sizeof(int)) ){

    //update the centroids
    updateCentroids(curr_assignment, points, cluster_centers, dim, N, nclusters);
    
    // for (int i = 0; i < nclusters; i++)
    //   printvec(&cluster_centers[i*dim], dim);

    memcpy(prev_assignment, curr_assignment, N*sizeof(int));

    clusterAssignment(points, curr_assignment, cluster_centers, dim, N, nclusters);
    epoch++;
    // printvec2(curr_assignment, N);
    cout<<"Epoch "<<epoch<<endl;
  }
  return epoch;
}

int main(int argc, char* argv[]){
  // vector<vector<double> > array;
  vector<double> finalarray;
  // vector<double> row;                /* vector for row of values */
  string line;

  int dim = 0;
  long N = 0;

  

  if (argc == 1){
    //random generate
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
        if (dim > 0)
          finalarray.insert(finalarray.end(), row.begin(), row.end());
    }

    N = finalarray.size();

    N = N/dim;
    
    double* finaldata = finalarray.data();

    /*
    cout << "complete array\n\n";
    for (int i = 0; i < N ; i++) {           
        for (int j = 0; j < dim; j++)
        cout << finaldata[i*dim + j] << "  ";
      cout<<endl;
    }*/


    // int nclusters = sscanf(argv[2]);
    int nclusters = 3;

    Timer tt;
    tt.tic();

    int *curr_assignment = (int *)calloc(N, sizeof(int));
    int epoch = trainKMeans(100, finaldata, N, dim, nclusters, curr_assignment);
    double elapsed = tt.toc();
    printf("Time elapsed  per epoch is %f seconds.\n", elapsed/epoch);
    
    // free(finaldata);
    f.close();
  }
  return 0;
}