compfile = input.huff
outfile = input.dec
file = input_1.txt 
threads = 10

pthreads: encode_pthreads decode_pthreads

omp: encode_openmp decode_openmp

encode_pthreads: encode_pthreads.c
	gcc -Wall  -o encode_pthreads.o encode_pthreads.c -lm -lpthread

decode_pthreads: decode_pthreads.c
	gcc -o decode_pthreads.o decode_pthreads.c -lm -lpthread

encode_openmp: encode_openmp.c
	gcc -Wall -o encode_openmp.o encode_openmp.c -lm -fopenmp

decode_openmp: decode_openmp.c
	gcc -o decode_openmp.o decode_openmp.c -lm -fopenmp

clean:
	rm *.o *.dec *.huff

pthreads_test: encode_pthreads decode_pthreads
	./encode_pthreads.o $(file) $(compfile) $(threads)
	./decode_pthreads.o $(compfile) $(outfile)
	diff $(file) $(outfile)

omp_test: encode_openmp decode_openmp
	./encode_openmp.o $(file) $(compfile)
	./decode_openmp.o $(compfile) $(outfile)
	diff $(file) $(outfile)