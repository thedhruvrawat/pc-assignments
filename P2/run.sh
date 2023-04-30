#!/bin/bash

file="input_1.txt" # Set the file to be tested here
compfile="input.huff"
outfile="input.dec"

i=1
while [ $i -le 10 ]  # continue while i is less than or equal to 10
do
  gcc -Wall -o encode_openmp.o encode_openmp.c -lm -fopenmp -D OMP_NUM_THREADS=$i
  gcc -o decode_openmp.o decode_openmp.c -lm -fopenmp -D OMP_NUM_THREADS=$i
  echo "Iteration $i"
  echo "Encoding and decoding with $i threads"
    compfile="input_${i}.huff"
    outfile="input_${i}.dec"
    ./encode_openmp.o $file $compfile $i
    ./decode_openmp.o $compfile $outfile $i
    diff $file $outfile
  i=$((i+1))  # increment i by 1
done

make clean