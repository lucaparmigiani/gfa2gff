SRC=main.cpp
all: ${SRC} 
	g++ -std=c++17 ${SRC} -fopenmp -lz -O2 -o gfa2gff
