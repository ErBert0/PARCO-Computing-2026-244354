#!/bin/bash


#Here you can modify the Compiler and the Compilation Flags

#module load gcc91

g++ -std=c++11 -O3 -o Sequential.out ../src/Sequential.cpp
g++ -std=c++11 -O3 -fopenmp -o Parallel.out ../src/Parallel.cpp

echo "Compilation Terminated."



