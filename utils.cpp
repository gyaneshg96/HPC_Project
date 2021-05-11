#include<iostream>
#include<cmath>
#include<string>

using namespace std;

double distance(double *x, double*y, int dim){
  //for now only euclidean distance
  double dist = 0;
  for (int i = 0; i < dim; i++){
    dist += (*x - *y)*(*x - *y);
    x++;
    y++;
  }
  return sqrt(dist);
}

//can use BLAS for this
void divide(double *x, double y, int dim){
  for (int i = 0; i < dim;i++){
    *x /= y;
    x++;
  }
}

void sum(double *x, double*y, int dim){
  for (int i = 0; i < dim;i++){
    *x += *y;
    x++;
    y++;
  }
}

void printvec(double *x, int dim){
  for (int i = 0; i <dim;i++){
    cout<<" "<<*x;
    x++;
  }
  cout<<endl;
}
void printvec2(int *x, int dim){
  for (int i = 0; i <dim;i++){
    cout<<" "<<*x;
    x++;
  }
  cout<<endl;
}