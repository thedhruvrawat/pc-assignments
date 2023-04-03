CC = mpicc
CFLAGS = -Wall -O2

all: blelloch hillis

blelloch: blelloch.c
	$(CC) $(CFLAGS) -o blelloch blelloch.c -lm

hillis: hillis.c
	$(CC) $(CFLAGS) -o hillis hillis.c -lm

clean:
	rm -f blelloch hillis *.out *.o *.txt