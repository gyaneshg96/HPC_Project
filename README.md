## Parallel KMeans with OpenMP, MPI, BLAS

This is a simple parallel KMeans implementation using multithreading, parallel processing using MPI
and high performance Linear Algebra library OpenBLAS. We use 3 different initialization, forgy, random, kmeans++.
We also use parallel file IO for improved performance.

### Requirements
 1. OpenMPI (for now we use the open source)
 2. OpenBLAS

###

### There are 3 files
 1. Simple serial KMeans
 2. Parallel KMeans with MPI/OpenMP
 3. Parallel KMeans with BLAS for vector operations

### There are 2 modes
 1. Mode 0 : specify folder containing csv files of the form 0.csv, 1.csv. Divide into chunks the file, for better preprocessing
 2. Mode 1 : specify the csv file 

### Usage

First, build the executables:

```make```

To run serial example, run

``` ./kmeans_serial <mode> <dataset> <nclusters> <initialization> ```

To run parallel run

```mpirun -np <nprocesses> ./kmeans_mpi <mode> <dataset> <nclusters> <initialization> ```


E.g.  

``` mpirun -np 4 ./kmeans_mpi 0 cifar/ 10 kmeans++ ```


