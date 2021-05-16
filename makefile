all:
	mpic++ -std=c++11 -O3 -march=native -fopenmp kmeansMPI.cpp -o kmeansMPI
