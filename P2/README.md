## Huffman Compression
This folder contains `pthreads` and `OpenMP` based implementations of Huffman encoding and decoding algorithms.

To run preset tests
```bash
bash run.sh
```

For individual results on specific number of processors (p)
```bash
make clean
# For pthreads
make pthreads
./encode_pthreads.o <input file> <compressed file> <no. of threads>
./decode_pthreads.o <compressed file> <output file>
diff <input file> <output file>
# For OpenMP
make omp
./encode_openmp.o <input file> <compressed file> <no. of threads>
./decode_openmp.o <compressed file> <output file> <no. of threads>
diff <input file> <output file>
```
