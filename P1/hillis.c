#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

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

int main(int argc, char** argv) {
    int data[MAX_SIZE];
    int element_count;
    read_input("input.txt", data, &element_count);
    int rank, size;
    int* local_data;
    int* local_y;
    int *prefix_sums = (int *)malloc(element_count * sizeof(int));

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);

    // if(rank==0) {   
    //     // Print the numbers to verify they were read correctly
    //     printf("Input:\t");
    //     for (int i = 0; i < element_count; i++)
    //     {
    //         printf("%d\t", data[i]);
    //     }
    //     printf("\n");
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

    // Allocate memory for local input array and output array
    local_data = (int*) malloc(sendcounts[rank] * sizeof(int));
    local_y = (int*) malloc(sendcounts[rank] * sizeof(int));
    
    MPI_Scatterv(data, sendcounts, displs, MPI_INT, local_data, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
    
    int num_steps = log2(element_count);
    int prefix_sum;

    double start_time = MPI_Wtime();

    // Perform Hillis-Steele parallel prefix sum
    for (int step = 0; step < num_steps; step++) {
        int power = pow(2, step);
        for (int k = 0; k < sendcounts[rank]; k++) {
            if (k >= power) {
                local_y[k] = local_data[k - power] + local_data[k];
            } else {
                local_y[k] = local_data[k];
            }
        }
        // Swap local_data and local_y
        int* temp = local_data;
        local_data = local_y;
        local_y = temp;
    }

    
    // Receive the prefix sum from the previous processes
    prefix_sum = 0;
    
    if(rank==0 && size>1) {
        MPI_Send(&local_data[sendcounts[rank] - 1], 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    if(rank>0 && rank!=size-1) {
        int t;
        MPI_Recv(&t, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        prefix_sum += t;
        for (int i = 0; i < sendcounts[rank]; i++) {
            local_data[i] += prefix_sum;
        }
        MPI_Send(&local_data[sendcounts[rank] - 1], 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    } 
    if(rank>0 && rank==size-1) {
        int t;
        MPI_Recv(&t, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        prefix_sum += t;
        for (int i = 0; i < sendcounts[rank]; i++) {
            local_data[i] += prefix_sum;
        }
    }

    double finish_time = MPI_Wtime(); 

    MPI_Gatherv(local_data, sendcounts[rank], MPI_INT, prefix_sums, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    // Print the prefix sums
    if (rank == 0) {
        // printf("Output:\t");
        // for (int i = 0; i < element_count; i++) {
        //     printf("%d\t", prefix_sums[i]);
        // }
        // printf("\n");
        printf("Parallel Time: %f\n", 1000 * (finish_time - start_time));
        free(prefix_sums);
    }

    // Clean up
    free(local_data);
    free(local_y);
    free(displs);
    free(sendcounts);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
