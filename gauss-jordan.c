#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "gauss-jordan.h"
#include "column.h"
#include "mpi_wrappers.h"

extern gauss_jordan gj;

// Forward declarations
void wait_all_wrapper(MPI_Request* requests, int size);
int is_my_column(int col_idx, int my_rank);
int process_of_column(int col_idx);
int number_of_cols(int rank);
int map_global_to_local(int k);
int number_of_recipients(int my_rank);
MPI_Request* send_column(int k, int max_idx, column_t* my_cols, int my_rank);
column_t receive_column(int k);

/**
 * Initializes the program, MPI framework, and sends the columns to the 
 * appropriate MPI processes
 * @param augmented_n The length of the 2D array (B-vector inclusive)
 * @param augmented_m The 2D array of coefficients
 * @param groupDistribution
 * @return The array of columns of each process
 */
column_t* init(int augmented_n, float** augmented_m, int useGroupDistribution) {
    gj.dimension = augmented_n - 1;
    gj.use_group_distribution = (useGroupDistribution != 0) ? 1 : 0;
    // Init MPI
    MPI_Init(NULL, NULL);
    // Initialize global Gauss Jordan info struct
    MPI_Comm_size(MPI_COMM_WORLD, &(gj.proc_num));
    gj.group_number = gj.dimension / gj.proc_num;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    column_t* my_cols;
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
        my_cols = malloc((gj.group_number + 1) * sizeof(column_t));
        MPI_Request* requests = 
            malloc(gj.group_number * (gj.proc_num - 1) * sizeof(MPI_Request));
        int my_cols_idx = 0, requests_idx = 0;
        for (j = 0; j < augmented_n - 1; j++) {
            int others_rank = process_of_column(j);
            if (others_rank != my_rank) {
                MPI_Isend(columns[j], gj.dimension, MPI_FLOAT, others_rank, j,
                    MPI_COMM_WORLD, &(requests[requests_idx++]));
            }
            else 
                my_cols[my_cols_idx++] = create_column(j, columns[j]);
        }
        // Add b-vector
        my_cols[my_cols_idx] = create_column(j, columns[gj.dimension]);

        // Wait on sent data
        wait_all_wrapper(requests, --requests_idx);
        free(requests);
        requests = NULL;
        // Free all columns array
        free(columns);
        columns = NULL;
//        printf("Master - Sent all data\n");
    }
    else {
        // Allocate space for your columns and receive them (synchronously)
        my_cols = malloc(gj.group_number * sizeof(column_t));
        float* tmp_data = malloc(gj.dimension * sizeof(float));
        int j, my_cols_idx = 0;
        for (j = 0; j < gj.dimension; j++) {
            if (is_my_column(j, my_rank)) {
//                printf("About to receive col[%d] -- col_idx == %d\n", j, my_cols_idx);
                MPI_Request request;
                MPI_Irecv(tmp_data, gj.dimension, MPI_FLOAT, 0,
                    j, MPI_COMM_WORLD, &request);
                wait_wrapper(&request);
//                printf("Slave(%d) - Got column[%d]:\n", my_rank, j);
                my_cols[my_cols_idx++] = create_column(j, tmp_data);
//                print_column(my_cols[my_cols_idx - 1]);
                if (my_cols_idx == gj.group_number)
                    break;
            }
        }
        free(tmp_data);
        tmp_data = NULL;
    }
    // Allocate space for the dummy column, used for sending and receiving
    // a column (acts as a serialization array)
    gj.dummy_col = malloc((gj.dimension + 1) * sizeof(float));
    free(augmented_m);
    return my_cols;
}

void destroy(column_t** my_cols_ptr, int my_rank) {
    // Free my_cols array
    int i;
    int upper_dim = (my_rank == 0) ? gj.group_number + 1 : gj.group_number;
    for (i = 0; i < upper_dim; i++) {
        delete_column((*my_cols_ptr)[i]);
    }
    free(*my_cols_ptr);
    *my_cols_ptr = NULL;
    // Free dummy column
    free(gj.dummy_col);
    gj.dummy_col = NULL;
}

void gj_kgi_main_loop(column_t* my_cols, int my_rank) {
    int k;
    for (k = 0; k < gj.dimension; k++) {
        if (is_my_column(k, my_rank)) {
//            printf("I'm rank(%d), and [%d] is my column\n", my_rank, k);
            column_t col = my_cols[map_global_to_local(k)];
            int pivot_idx = pivot(col, k);
            MPI_Request* requests = send_column(k, pivot_idx, my_cols, my_rank);
            int j, upper_bound = number_of_cols(my_rank);
            for (j = map_global_to_local(k) + 1; j < upper_bound; j++) {
//                if (my_rank == 0) 
//                    printf("%d-Proc[%d] modifying column[%d]\n", k, my_rank, j);
                modify(my_cols[j], col, pivot_idx, k);
            }
            int recipients_number = (my_rank == gj.proc_num - 1) ? 
                1 : number_of_recipients(my_rank);
//            printf("\t%d-Proc[%d] WILL WAIT\n", k, my_rank);
            wait_all_wrapper(requests, recipients_number);
//            printf("__\t%d-Proc[%d] finished waiting\n", k, my_rank);
            free(requests);
        }
        else {
            // Receive column and extract pivot index
            if ((my_rank == 0) || 
                    ( !(gj.use_group_distribution && (my_rank < process_of_column(k))) )) {
//                printf("\t%d-Proc[%d] WILL RECEIVE col\n", k, my_rank);
                column_t col = receive_column(k);
//                printf("\t\t%d-Proc[%d] RECEIVED col\n", k, my_rank);
                int pivot_idx = (int) gj.dummy_col[0];
                /* 
                 * Which local columns to modify?
                 *      On k-group distribution, start from 0. As all the previous,
                 *          columns, belong to my_rank - 1 processes
                 *      On circular distribution, start from the conversion of 
                 *          k to the local_idx format + 1. Although, if the 
                 *          process_of_col(k) isn't the last one, the processes that
                 *          receive the k-column and have rank > my_rank must modify
                 *          the local_idx local column, as well
                 */
                if (gj.use_group_distribution) {
                    if (my_rank == 0) {
//                        printf("\t%d--Proc[0] modifying column[%d]\n", k, gj.group_number);
                        modify(my_cols[gj.group_number], col, pivot_idx, k);
                    }
                    else {
                        int j;
                        for (j = 0; j < gj.group_number; j++) {
//                            if (my_rank == 0) 
//                                printf("\t%d--Proc[%d] modifying column[%d]\n", k, my_rank, j);
                            modify(my_cols[j], col, pivot_idx, k);
                        }
                    }
                }
                else {
                    int local_idx = map_global_to_local(k);
                    // Not the latest column, start modifying from the local_idx
                    int j, start_idx = (process_of_column(k) != gj.proc_num - 1) ? 
                            local_idx : local_idx + 1;
                    int upper_bound = (my_rank == 0) ? 
                        gj.group_number + 1 : gj.group_number;
                    for (j = start_idx; j < upper_bound; j++)
                        modify(my_cols[j], col, pivot_idx, k);
                }
                delete_column(col);
            }
        }
    }
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
    if (gj.use_group_distribution)
        return col_idx / gj.group_number;
    else
        return col_idx % gj.proc_num;
}

/**
 * Returns the number of columns that a process has
 * (group number := dimension/process_number)
 * Rank_0 has group_number + 1, all others group_number
 * @param rank The rank of the process whose number of columns are needed
 * @return The number of columns of the rank-th process
 */
int number_of_cols(int rank) {
    return (rank == 0) ? gj.group_number + 1 : gj.group_number;
}

/**
 * Maps the global idx (k - step) to the idx of the local array
 * @param k The k-step idx
 * @return The local idx
 */
int map_global_to_local(int k) {
    return ((gj.use_group_distribution) ? k % gj.group_number : k / gj.proc_num);
}

/**
 * Returns the number of recipients depending on the distribution
 * @param my_rank The rank of the current MPI Process
 * @return The number of recipients
 */
int number_of_recipients(int my_rank) {
    return (gj.use_group_distribution) ? gj.proc_num - my_rank - 1 : gj.proc_num - 1;
}

/**
 * Send the k-th column to all the > k column owners
 * @param k The current step of the method
 * @param max_idx The pivot element index
 * @param my_cols The columns array of the current MPI process
 * @param my_rank The rank of the current MPI Process
 * @return Array of MPI_Requests, that was populated on Isend
 */
MPI_Request* send_column(int k, int max_idx, column_t* my_cols, int my_rank) {
    // Get the k-th column
    int local_idx = map_global_to_local(k);
    // Serialize the struct to the local array
    gj.dummy_col[0] = max_idx;
    memcpy(&(gj.dummy_col[1]), my_cols[local_idx]->data, gj.dimension * sizeof(float));
    // Send it to all the > k owners, but gather indices to my own
    // (that are > k)
    // No need to send columns to others except for master process (has b-vector) 
    // if in group-distribution and you're the latest process
    if (gj.use_group_distribution && (my_rank == gj.proc_num - 1)) {
        MPI_Request *requests = malloc(sizeof(MPI_Request));
        MPI_Isend(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, 0, k, 
                    MPI_COMM_WORLD, &(requests[0]));
        return requests;
    }
    MPI_Request *requests = malloc(
            number_of_recipients(my_rank) * sizeof(MPI_Request));
    int rank, requests_idx = 0;
    if (gj.use_group_distribution) {
        for (rank = my_rank + 1; rank < gj.proc_num; rank++) {
//            printf("Proc[%d] sending to proc[%d]\n", my_rank, rank);
            MPI_Isend(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, rank, k, 
                    MPI_COMM_WORLD, &(requests[requests_idx++]));
        }
        // Send to master process as well
        MPI_Isend(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, 0, k, 
                    MPI_COMM_WORLD, &(requests[requests_idx++]));
    }
    else {
        for (rank = 0; rank < gj.proc_num; rank++) {
            if (rank != my_rank) 
                MPI_Isend(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, rank, k, 
                    MPI_COMM_WORLD, &(requests[requests_idx++]));
        }
    }
//    if (k == 1) {
//        printf("[k==1], sent to %d processes\n", requests_idx);
//    }
    return requests;
}

/**
 * Receives the column k and saves it to the dummy_col field
 * @param k The current step of the Gauss Jordan method
 * @return The newly arrived column
 */
column_t receive_column(int k) {
    MPI_Request request;
    MPI_Irecv(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, process_of_column(k), 
            k, MPI_COMM_WORLD, &request);
    wait_wrapper(&request);
    return create_column(k, &(gj.dummy_col[1]));
}

/**
 * Randomly creates a [n x n] matrix
 * @param dimension The 1D of the matrix
 * @return The matrix
 */
float** create_square_matrix(int dimension) {
    int i,j;
    float** matrix = malloc(dimension * sizeof(float*));
    for (i = 0; i < dimension; i++) {
        matrix[i] = malloc(dimension * sizeof(float));
        for (j = 0; j < dimension; j++) {
            matrix[i][j] = rand();
        }
    }
    return matrix;
}

float* create_random_array(int dimension) {
    float* array = malloc(dimension * sizeof(float));
    int i;
    for (i = 0; i < dimension; i++) {
        array[i] = rand();
    }
    return array;
}