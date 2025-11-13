#!/bin/bash

module load gcc91

g++ -std=c++11 -march=native -O3 -o Sequential.out Sequential.cpp

g++ -std=c++11 -march=native -O3 -fopenmp -o Parallel.out Parallel.cpp

echo "Compilazione terminata. Eseguibili creati:"
ls -lh *.out
