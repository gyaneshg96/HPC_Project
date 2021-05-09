#include<iostream>
#include<string>
#include"utils.cpp"
#include<vector>
#include<fstream>
#include<cstring>
#include<sstream>

using namespace std;



double ** initialize(double **points, long N, int dim, int nclusters, string method){
  double **cluster_centers = (double **)malloc(nclusters*sizeof(double *));
  for (int i = 0; i < nclusters; i++){
    cluster_centers[i] = (double *)calloc(dim, sizeof(double));
  }

  if (method == "forgy"){
    for (int i = 0; i < nclusters;i++){
      int val = rand() % N;
      memcpy(cluster_centers[i], points[val], dim*sizeof(double));
      //which clusters
  }
  }
  if (method == "random"){
    int *cluster_nums = (int *)calloc(nclusters, sizeof(int));
    for (long i = 0; i< N;i++){
      int cluster = rand() % nclusters;
      //add to cluster centers
      // cluster_centers[cluster] += points[i];
      sum(cluster_centers[cluster], points[i], dim);
      cluster_nums[cluster] += 1;
    }
    for (int i = 0; i < nclusters; i++)
      // cluster_centers[i] /= cluster_nums[i];
      divide(cluster_centers[i], cluster_nums[i], dim);
    
    free(cluster_nums);
  }
  /*
  if (method == "kmeans++"){
    int val = rand() % N;
    cluster_centers[0] = points[val];
    for (int i = 1; i < nclusters; i++){
      double mindist = 100000;
      for (int j = 0; j < N;i++){
        mindist = min(mindist, distance(points[j], cluster_centers[i], dim));
    }
  }
  */
  return cluster_centers;
}

void clusterAssignment(double **points, int* cluster_assignment, double **cluster_centers, int dim, long N, int nclusters){

  for (long point = 0; point < N; point++){

    double min_distance = 100000;
    double min_center = 0;
    double *currpoint = points[point];    

    for (int center = 0; center < nclusters; center++){
      double dist = distance(cluster_centers[center], points[point], dim);
      if (dist < min_distance){
        min_distance = dist;
        min_center = center; 
      }
    }
    // cout<<min_center<<endl;
    cluster_assignment[point] = min_center;
  }
}
void updateCentroids(int *cluster_assignment, double **points, double **cluster_centers, int dim, long N, int nclusters){
  
  for (int i = 0; i < nclusters; i++){
    for(int j = 0; j < dim; j++)
      cluster_centers[i][j] = 0; 
  }

  int *cluster_nums = (int *)calloc(nclusters, sizeof(int));
  for (int point = 0; point < N; point++){
    sum(cluster_centers[cluster_assignment[point]], points[point], dim);
    cluster_nums[cluster_assignment[point]] += 1;    
  }
  // printvec(cluster_centers[0], dim);

  for (int i = 0; i < nclusters; i++){
    divide(cluster_centers[i], cluster_nums[i], dim);
  }

  free(cluster_nums);
}

//if prev assignment is the same as current
int compareAssignment(int *prevassignment, int *currassignment, long N){
  for (long i = 0; i < N;i++ ){
    if (prevassignment[i] != currassignment[i])
      return 1;
  }
  return 0;
}
void trainKMeans(int epochs, double **points, long N, int dim, int nclusters){  
  int epoch = 0;
  double **cluster_centers = initialize(points, N, dim, nclusters, "random");

  // for (int i = 0; i < nclusters; i++)
    // printvec(cluster_centers[i], dim);

  cout<<"Init"<<endl;

  int *curr_assignment = (int *)calloc(N, sizeof(int));
  clusterAssignment(points, curr_assignment, cluster_centers, dim, N, nclusters);

  cout<<"Cluster Assigned"<<endl;
  //all zero
  // printvec2(curr_assignment, N);
  int *prev_assignment = (int *)calloc(N, sizeof(int));

  // while(compareAssignment(prev_assignment, curr_assignment, N) && epoch < epochs){
  while(memcmp(prev_assignment, curr_assignment, N*sizeof(int)) && epoch < epochs){

    //update the centroids
    updateCentroids(curr_assignment, points, cluster_centers, dim, N, nclusters);
    
    for (int i = 0; i < nclusters; i++)
      printvec(cluster_centers[i], dim);

    memcpy(prev_assignment, curr_assignment, N*sizeof(int));

    clusterAssignment(points, curr_assignment, cluster_centers, dim, N, nclusters);
    epoch++;
    // printvec2(curr_assignment, N);
    cout<<"Epoch "<<epoch<<endl;
  }
}

int main(int argc, char* argv[]){
  vector<vector<double> > array;
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
        vector<double> row;                /* vector for row of values */
        stringstream s (line);        /* stringstream to parse csv */
        while (getline(s, val, ',')){
            try {
            double dval = stod(val); 
            row.push_back(stod(val));
            }
            catch (exception& e){
              //do nothing
            }
            }
        dim = row.size();
        if (dim > 0)
          array.push_back(row);         /* add row to array */
    }
    N = array.size();
    
    double **finalarray = (double **)malloc(N*sizeof(double *));
    for (int i = 0; i< N;i++){
      finalarray[i] = array[i].data();
    }


    /*
    cout << "complete array\n\n";
    for (int i = 0; i < N ; i++) {           
        for (int j = 0; j < dim; j++)
        cout << finalarray[i][j] << "  ";
    }*/

    cout<<N<<" "<<dim<<endl;
    // int nclusters = sscanf(argv[2]);
    int nclusters = 3;
    trainKMeans(100, finalarray, N, dim, nclusters);
    f.close();
  }
  return 0;
}