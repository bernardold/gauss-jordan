#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "gauss-jordan.h"
#include "column.h"
#include "mpi_wrappers.h"

extern gauss_jordan gj;

// Forward declarations
void wait_all_wrapper(MPI_Request* requests, int size);
int is_my_column(int col_idx, int my_rank);
int process_of_column(int col_idx);

/**
 * Initializes the program, MPI framework, and sends the columns to the 
 * appropriate MPI processes
 * @param augmented_n The length of the 2D array (B-vector inclusive)
 * @param augmented_m The 2D array of coefficients
 * @param groupDistribution
 * @return The array of columns of each process
 */
column** init(int augmented_n, float** augmented_m, int groupDistribution) {
    gj.dimension = augmented_n - 1;
    gj.groupDistribution = groupDistribution;
    // Init MPI
    MPI_Init(NULL, NULL);
    // Initialize global Gauss Jordan info struct
    MPI_Comm_size(MPI_COMM_WORLD, &(gj.proc_num));
    gj.groupNumber = gj.dimension / gj.proc_num;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    column** my_cols;
    // Master populates an all-columns array
    // then sends the appropriate columns to their "owners"
    // Others, receive their columns and add them to their column array
    if (my_rank == 0) {
        // Convert array to columns
        float** columns = malloc(augmented_n * sizeof(float*));
        int j;
        for (j = 0; j < augmented_n; j++) {
            columns[j] = malloc((augmented_n - 1) * sizeof(float));
            int i;
            for (i = 0; i < augmented_n - 1; i++) {
                columns[j][i] = augmented_m[i][j];
            }
        }
        
        // Send to others and populate yours
        // (master handles the b-vector as well)
        my_cols = malloc((gj.groupNumber + 1) * sizeof(column*));
        MPI_Request* requests = 
            malloc(gj.groupNumber * (gj.proc_num - 1) * sizeof(MPI_Request));
        int my_cols_idx = 0, requests_idx = 0;
        for (j = 0; j < augmented_n - 1; j++) {
            int others_rank = process_of_column(j);
            if (others_rank != my_rank) {
                MPI_Isend(columns[j], gj.dimension, MPI_FLOAT, others_rank, j,
                    MPI_COMM_WORLD, &(requests[requests_idx++]));
            }
            else 
                my_cols[my_cols_idx] = create_column(j, columns[j]);
        }
        // Add b-vector
        my_cols[my_cols_idx] = create_column(j, columns[gj.dimension]);

        // Wait on sent data
        wait_all_wrapper(requests, --requests_idx);
        // Free all columns array
        printf("Master - Sent all data\n");
    }
    else {
        // Allocate space for your columns and receive them (synchronously)
        my_cols = malloc(gj.groupNumber * sizeof(column));
        float* tmp_data = malloc(gj.dimension * sizeof(float));
        int j, my_cols_idx = 0;
        for (j = 0; j < gj.dimension; j++) {
            if (is_my_column(j, my_rank)) {
                printf("About to receive col[%d] -- col_idx == %d\n", j, my_cols_idx);
                MPI_Request request;
                MPI_Irecv(tmp_data, gj.dimension, MPI_FLOAT, 0,
                    j, MPI_COMM_WORLD, &request);
                wait_wrapper(&request);
                printf("Slave(%d) - Got column[%d]:\n", my_rank, j);
                my_cols[my_cols_idx] = create_column(j, tmp_data);
                print_column(my_cols[my_cols_idx]);
                if (my_cols_idx == gj.groupNumber)
                    break;
            }
        }
        free(tmp_data);
        tmp_data = NULL;
    }
    return my_cols;
}

/**
 * Checks if the column given belongs to the MPI process passed
 * @param col_idx The column idx
 * @param my_rank The MPI Process to be checked against owning the column
 * @return On true, 1 is returned, otherwise 0
 */
int is_my_column(int col_idx, int my_rank) {
    return process_of_column(col_idx) == my_rank;
}

/**
 * Returns the ID (rank) of the MPI Process that owns the given column id
 * @param col_idx The column id (from 0..n)
 * @return The rank of the MPI Process
 */
int process_of_column(int col_idx) {
    if (gj.groupDistribution) {
        return col_idx / gj.groupNumber;
    }
    else {
        return col_idx % gj.proc_num;
    }
}