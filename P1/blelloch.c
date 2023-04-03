#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>

#define MAX_SIZE 2000001
#define INPUT_FILE "input.txt"

void read_input(char* file_name, int* data, int* n) {
    FILE* fp = fopen(file_name, "r");
    if (fp == NULL) {
        printf("Error opening file %s.\n", file_name);
        return;
    }
    if (fscanf(fp, "%d", n) != 1) {
        printf("Error reading input!\n");
        exit(1);
    }
    for (int i = 0; i < *n; i++) {
        if (fscanf(fp, "%d", &data[i]) != 1) {
            printf("Error reading input!\n");
            exit(1);
        }
    }
    fclose(fp);
}

int getSum(int* input, size_t count) {
    int result = 0;
    for (int i = 0; i<count; i++) {
        result += input[i];
    }
    // printf("*%d*", result);
    return result;
}

void find_sum(int rank, int size, int *input, size_t elements_per_process, int* total) {
    MPI_Status status;

    int sum = getSum(input, elements_per_process);

    bool isAlive = true;

    for (int level = 0; level < (int)log2(size); level++) {
        if(isAlive) {
            int pos = rank / (int)pow(2, level);
            if(pos%2==0) {
                int temp = 0;
                int getFrom = rank + (int)pow(2, level);
                MPI_Recv(&temp, 1, MPI_INT, getFrom, 0, MPI_COMM_WORLD, &status);
                sum+=temp;
            } else {
                int sendTo = rank - (int)pow(2, level);
                MPI_Send(&sum, 1, MPI_INT, sendTo, 0, MPI_COMM_WORLD);
                isAlive = false;
            } 
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    *total = sum;
}

void find_prefix_sum(int rank, int size, int *input, size_t elements_per_process, int sum) {
    int prefixSum = 0;
    MPI_Status status;
    if(rank == 0) {
        prefixSum = sum;
    }

    for (int level = (int)log2(size) - 1; level>=0; level--) {
        // Only check processes on this level
        if(level == 0 || rank%(int)pow(2,level) == 0) {
            int pos = rank / (int)pow(2, level);
            if(pos%2==1) {
                int parent = rank - (int)pow(2, level);
                MPI_Recv(&prefixSum, 1, MPI_INT, parent, 0, MPI_COMM_WORLD, &status);
                MPI_Send(&sum, 1, MPI_INT, parent, 0, MPI_COMM_WORLD);
            } else {
                int temp = 0;
                int rightChild = rank + (int)pow(2, level);
                MPI_Send(&prefixSum, 1, MPI_INT, rightChild, 0, MPI_COMM_WORLD);
                MPI_Recv(&temp, 1, MPI_INT, rightChild, 0, MPI_COMM_WORLD, &status);
                prefixSum-=temp;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int next_sum = input[elements_per_process - 1];
    // printf("%d ", next_sum);
    input[elements_per_process - 1] = prefixSum;
    for (int j = elements_per_process - 2; j >= 0; j--) {
        int temp = input[j];
        input[j] = input[j + 1] - next_sum;
        next_sum = temp;
    }
}

int main(int argc, char *argv[]) {
    int data[MAX_SIZE];
    int element_count;
    read_input("input.txt", data, &element_count);
    
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);

    // if(rank==0) {   
        // Print the numbers to verify they were read correctly
        // printf("Input:\t");
        // for (int i = 0; i < element_count; i++)
        // {
        //     printf("%d\t", data[i]);
        // }
        // printf("\n");
        // start_time_seq = MPI_Wtime();
        // int* result_seq = NULL;
        // result_seq = malloc(element_count * sizeof(int));
        // prefixSum_sequential(data, result_seq, element_count);
        // finish_time_seq = MPI_Wtime(); 
    // }
    MPI_Barrier(MPI_COMM_WORLD);
    size_t elements_per_process = element_count / size;
    int remainder = element_count % size;
    int* sendcounts = malloc(size * sizeof(int));
    int* displs = malloc(size * sizeof(int));
    for (int i = 0; i < size; i++) {
        sendcounts[i] = elements_per_process;
        if (i < remainder) {
            sendcounts[i]++;
        }
        displs[i] = i * elements_per_process + ((i < remainder) ? i : remainder);
    }
    int* local_data = malloc(sendcounts[rank] * sizeof(int));

    MPI_Scatterv(data, sendcounts, displs, MPI_INT, local_data, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
    int sum;

    double start_time = MPI_Wtime();

    find_sum(rank, size, local_data, sendcounts[rank], &sum);
    find_prefix_sum(rank, size, local_data, sendcounts[rank], sum);

    double finish_time = MPI_Wtime(); 

    int* result = NULL;
    if (rank == 0) {
        result = malloc(element_count * sizeof(int));
    }
    MPI_Gatherv(local_data, sendcounts[rank], MPI_INT, result, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // printf("Output:\t");
        // for (int i = 0; i < element_count; i++)
        // {
        //     printf("%d\t", result[i]);
        // }
        // printf("\n");
        printf("Parallel Time: %f\n", 1000 * (finish_time - start_time));
        free(result);
    }
    free(local_data);
    free(displs);
    free(sendcounts);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}