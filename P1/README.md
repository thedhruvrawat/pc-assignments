## Prefix Sum Algorithms
This folder contains `MPI` based implementations of Blelloch's Prefix Sum algorithm and Hillis-Steele Prefix Sum algorithm.

To run preset tests
```bash
bash run.sh
```

For individual results on specific number of processors (p)
```bash
make clean
# Generate n random numbers
python3 generate.py <n>
# For Hillis-Steele Algorithm
make hillis
mpirun -np <p> ./hillis
# For Blelloch's Algorithm
make blelloch
mpirun -np <p> ./blelloch
```
