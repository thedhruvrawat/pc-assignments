#!/bin/bash

make clean

make all
echo "-------------------------------------------------------------------"
echo "Testing for n=1000"
# Generate 1000 random integers using Python script
python3 generate.py 1000

echo "p=1"
mpirun -np 1 ./blelloch
echo "p=2"
mpirun -np 2 ./blelloch
echo "p=4"
mpirun -np 4 ./blelloch
echo "p=8"
mpirun -np 8 ./blelloch
echo "p=16"
mpirun -np 16 ./blelloch
echo "p=32"
mpirun -np 32 ./blelloch
echo "-------------------------------------------------------------------"
echo "Testing for n=10000"
# Generate 1000 random integers using Python script
python3 generate.py 10000

echo "p=1"
mpirun -np 1 ./blelloch
echo "p=2"
mpirun -np 2 ./blelloch
echo "p=4"
mpirun -np 4 ./blelloch
echo "p=8"
mpirun -np 8 ./blelloch
echo "p=16"
mpirun -np 16 ./blelloch
echo "p=32"
mpirun -np 32 ./blelloch
echo "-------------------------------------------------------------------"
echo "Testing for n=100000"
# Generate 1000 random integers using Python script
python3 generate.py 100000

echo "p=1"
mpirun -np 1 ./blelloch
echo "p=2"
mpirun -np 2 ./blelloch
echo "p=4"
mpirun -np 4 ./blelloch
echo "p=8"
mpirun -np 8 ./blelloch
echo "p=16"
mpirun -np 16 ./blelloch
echo "p=32"
mpirun -np 32 ./blelloch
echo "-------------------------------------------------------------------"
echo "Testing for n=500000"
# Generate 1000 random integers using Python script
python3 generate.py 500000

echo "p=1"
mpirun -np 1 ./blelloch
echo "p=2"
mpirun -np 2 ./blelloch
echo "p=4"
mpirun -np 4 ./blelloch
echo "p=8"
mpirun -np 8 ./blelloch
echo "p=16"
mpirun -np 16 ./blelloch
echo "p=32"
mpirun -np 32 ./blelloch

echo "-------------------------------------------------------------------"
echo "Testing for n=1000000"
# Generate 1000 random integers using Python script
python3 generate.py 1000000

echo "p=1"
mpirun -np 1 ./blelloch
echo "p=2"
mpirun -np 2 ./blelloch
echo "p=4"
mpirun -np 4 ./blelloch
echo "p=8"
mpirun -np 8 ./blelloch
echo "p=16"
mpirun -np 16 ./blelloch
echo "p=32"
mpirun -np 32 ./blelloch

echo "-------------------------------------------------------------------"
echo "Testing for n=2000000"
# Generate 1000 random integers using Python script
python3 generate.py 2000000

echo "p=1"
mpirun -np 1 ./blelloch
echo "p=2"
mpirun -np 2 ./blelloch
echo "p=4"
mpirun -np 4 ./blelloch
echo "p=8"
mpirun -np 8 ./blelloch
echo "p=16"
mpirun -np 16 ./blelloch
echo "p=32"
mpirun -np 32 ./blelloch
echo "-------------------------------------------------------------------"
echo "Testing for n=1000"
# Generate 1000 random integers using Python script
python3 generate.py 1000

echo "p=1"
mpirun -np 1 ./hillis
echo "p=2"
mpirun -np 2 ./hillis
echo "p=4"
mpirun -np 4 ./hillis
echo "p=8"
mpirun -np 8 ./hillis
echo "p=16"
mpirun -np 16 ./hillis
echo "p=32"
mpirun -np 32 ./hillis
echo "-------------------------------------------------------------------"
echo "Testing for n=10000"
# Generate 1000 random integers using Python script
python3 generate.py 10000

echo "p=1"
mpirun -np 1 ./hillis
echo "p=2"
mpirun -np 2 ./hillis
echo "p=4"
mpirun -np 4 ./hillis
echo "p=8"
mpirun -np 8 ./hillis
echo "p=16"
mpirun -np 16 ./hillis
echo "p=32"
mpirun -np 32 ./hillis
echo "-------------------------------------------------------------------"
echo "Testing for n=100000"
# Generate 1000 random integers using Python script
python3 generate.py 100000

echo "p=1"
mpirun -np 1 ./hillis
echo "p=2"
mpirun -np 2 ./hillis
echo "p=4"
mpirun -np 4 ./hillis
echo "p=8"
mpirun -np 8 ./hillis
echo "p=16"
mpirun -np 16 ./hillis
echo "p=32"
mpirun -np 32 ./hillis
echo "-------------------------------------------------------------------"
echo "Testing for n=500000"
# Generate 1000 random integers using Python script
python3 generate.py 500000

echo "p=1"
mpirun -np 1 ./hillis
echo "p=2"
mpirun -np 2 ./hillis
echo "p=4"
mpirun -np 4 ./hillis
echo "p=8"
mpirun -np 8 ./hillis
echo "p=16"
mpirun -np 16 ./hillis
echo "p=32"
mpirun -np 32 ./hillis

echo "-------------------------------------------------------------------"
echo "Testing for n=1000000"
# Generate 1000 random integers using Python script
python3 generate.py 1000000

echo "p=1"
mpirun -np 1 ./hillis
echo "p=2"
mpirun -np 2 ./hillis
echo "p=4"
mpirun -np 4 ./hillis
echo "p=8"
mpirun -np 8 ./hillis
echo "p=16"
mpirun -np 16 ./hillis
echo "p=32"
mpirun -np 32 ./hillis

echo "-------------------------------------------------------------------"
echo "Testing for n=2000000"
# Generate 1000 random integers using Python script
python3 generate.py 2000000

echo "p=1"
mpirun -np 1 ./hillis
echo "p=2"
mpirun -np 2 ./hillis
echo "p=4"
mpirun -np 4 ./hillis
echo "p=8"
mpirun -np 8 ./hillis
echo "p=16"
mpirun -np 16 ./hillis
echo "p=32"
mpirun -np 32 ./hillis
echo "-------------------------------------------------------------------"