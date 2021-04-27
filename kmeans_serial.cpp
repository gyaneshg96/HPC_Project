#include<iostream>
#include<string>
#include"utils.cpp"

using namespace std;

void initialize(double *points[], long N, int nclusters, string method){
  int *cluster = (int *)malloc(N*sizeof(int));
  double *cluster_centers[] = (int *)malloc(nclusters*sizeof(int));

  if (method == "forgy"){
    for (int i = 0; i < nclusters;i++){
      cluster_centers[i] = rand() % N;
      //which clusters
  }
  if (method == "random"){
    
  }
  if (method == "kmeans++")

}
void clusterAssignment(double *points[], double *cluster_centers[], long N, int clusters){
  int* cluster_assignment = (int *)malloc(N*sizeof(int));

  for (int point = 0; point < N; point++){
    double min_distance = INT16_MAX*1.0;
    double min_center = 0;
    for (int center = 0; center < nclusters; center++){
      double dist = distance(cluster_centers[center], points[points]);
      if (dist < min_distance){
        min_distance = dist;
        min_center = center; 
      }
    }
    cluster_assignment[point] = min_center;
  }
}
void findMean(int *cluster_assignment, double *points[], double *cluster_centers[], long N, int nclusters){
  double *cluster_means[] = {}
  for (int point = 0; point < N; point++){
    int 
  }  
}
void trainKMeans(int epochs, double error, ){
  int epochs = 1000
  while()
}