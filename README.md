# CS F422: Parallel Computing Assignment

## P1
Implement **Blelloch’s scan algorithm** and **Hillis and Steele’s algorithm** using `MPI`. You can take a list of numbers in a file `input.txt`.

1. Draw a task dependency graph for the parallel tasks. Identify opportunities for data paralelism, functional parallelism or pipelining. What is degree of concurrency?

2. Explain your design in translating the identified parallelism into `MPI`.

3. Compute speedup, efficiency, cost and isoefficiency metrics in terms of `n`, and `p` where `n` is number of data elements and `p` is number of processors.

4. Findout whether these algorithms are cost-optimal?

5. Evaluate the speedup and efficiency by running your program for a range of processes.

### Deliverables: 
- Design Document (`.pdf`). Must contain answers for 1-5. 
- Source code `blelloch.c` and `hillis.c`

## P2
Consider a text file (of atleast 1 MB in size) to be encoded using **Huffman codes**. Now consider parallel algorithms for encoding a given text file and decoding a given encoded file respectively.

1. Draw a task dependency graph for the parallel tasks. Identify opportunities for data paralelism, functional parallelism or pipelining. What is degree of concurrency?

2. Explain your design in translating the identified parallelism into `Pthreads`.

3. Compute speedup, efficiency, cost and isoefficiency metrics in terms of `n`, and `p` where `n` is number of data elements and `p` is number of processors.

4. Findout whether these algorithms are cost-optimal?

5. Evaluate the speedup and efficiency by running your program for a range of processes.

6. Using `Pthreads`, devise and implement algorithms for encoding and decoding.

7. Taking any open source serial implementation, apply `OpenMP` directives for encoding and decoding.

8. Compare speedup between 6 and 7

### Deliverables: 
- Design Document (`.pdf`) explaining your design. Must contain answers for 1-5 and 8. 
- Source code for 6: `encode_pthreads.c`, `decode_pthreads.c`
- Source code for 7: `encode_openmp.c`, `decode_openmp.c`

## P3
Consider a problem of printing all **Kaprekar numbers** upto given N taken as CLA. What is Kaprekar number? Consider an n-digit number k. Square it and add the right n digits to the left n or n-1 digits. If the resultant sum is k, then k is called a Kaprekar number. For example, 9 is a Kaprekar number since $9^2 = 81$, and $8 + 1 = 9$. $297^2 = 88209$ and $88 + 209 = 297$.

1. Draw a task dependency graph for the parallel tasks. Identify opportunities for data paralelism, functional parallelism or pipelining. What is degree of concurrency?

2. Explain your design in translating the identified parallelism into `CUDA`.

3. Compute speedup, efficiency, cost and isoefficiency metrics in terms of `n`, and `p` where `n` is number of data elements and `p` is number of processors.

4. Findout whether these algorithms are cost-optimal?

5. Evaluate the speedup and efficiency by running your program for a range of processes.

6. Implement your algorithm in `CUDA`

7. Taking any open source serial implementation, apply `OpenACC` directives for computing Kaprekar numbers.

8. Compare speedup between 6 and 7

### Deliverables: 
- Design Document (`.pdf`). Must contain answers for 1-5 and 8. 
- Source code: `kaprekar_cuda.cu`
- Source code: `kaprekar_openacc.c`