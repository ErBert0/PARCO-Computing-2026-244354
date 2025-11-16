#!/bin/bash


#This finds the root directory of the project
ROOT=$(cd "$(dirname "$0")/.." && pwd)

module load gcc91

#Here you can modify the Compiler and the Compilation Flags
COMPILER="g++-9.1.0"
FLAGS="-O3 -fopenmp"

$COMPILER $FLAGS -o "$ROOT/Sequential.out" "$ROOT/src/Sequential.cpp"
$COMPILER $FLAGS -o "$ROOT/Parallel.out" "$ROOT/src/Parallel.cpp"

echo "Compilation Terminated."



