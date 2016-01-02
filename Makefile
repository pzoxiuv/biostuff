all:
	g++ -Wall -O3 -std=c++11 -fopenmp main.cpp encoding.cpp entropy.cpp
