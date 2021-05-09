#include<iostream>
#include<cmath>
#include<string>

using namespace std;

double distance(double *x, double*y, int dim){
  //for now only euclidean distance
  double dist = 0;
  for (int i = 0; i < dim; i++){
    dist += (x[i] - y[i])*(x[i] - y[i]);
  }
  return sqrt(dist);
}

//can use BLAS for this
void divide(double *x, double y, int dim){
  for (int i = 0; i < dim;i++){
    x[i] /= y;
  }
}

void sum(double *x, double*y, int dim){
  for (int i = 0; i < dim;i++){
    x[i] += y[i];
  }
}

void printvec(double *x, int dim){
  for (int i = 0; i <dim;i++){
    cout<<" "<<x[i];
  }
  cout<<endl;
}
void printvec2(int *x, int dim){
  for (int i = 0; i <dim;i++){
    cout<<" "<<x[i];
  }
  cout<<endl;
}