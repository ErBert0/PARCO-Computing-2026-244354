#!/bin/bash


#This finds the root directory of the project
ROOT=$(cd "$(dirname "$0")/.." && pwd)

module load gcc91

#Here you can modify the Compiler and the Compilation Flags
COMPILER="g++-9.1.0"

$COMPILER -O3 -march=native -o "$ROOT/Sequential.out" "$ROOT/src/Sequential.cpp"
$COMPILER -O3 -march=native -fopenmp -o "$ROOT/Parallel.out" "$ROOT/src/Parallel.cpp"

echo "Compilation Terminated."