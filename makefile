all:
	mpic++ -std=c++11 -O3 -march=native -fopenmp -g kmeansMPI.cpp -o kmeansMPI
