all:
	mpic++ -std=c++11 -O3 -march-native -fopenmp -g kmeansMPI.cpp -o kmeansMPI
	g++ -std=c++11 -O3 kmeans_serial.cpp -o kmeans_serial
	mpic++ -std=c++11 -O3 -fopenmp  kmeans_mpi_blas.cpp -I /opt/OpenBLAS/include/ -L /opt/OpenBLAS/lib/ -lopenblas -o kmeans_mpi_blas
	mpic++ -std=c++11 -O3 -fopenmp  kmeans_mpi.cpp -o kmeans_mpi
clean:
	rm kmeans_serial kmeans_mpi kmeans_mpi_blas
