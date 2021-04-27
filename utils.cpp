#include<iostream>
#include<math>
#include<string>

double distance(double *x, double*y, int N){
  //for now only euclidean distance
  double dist = 0;
  for (int i = 0; i < N; i++){
    dist += x[i]*y[i];
  }
  return sqrt(dist);
}