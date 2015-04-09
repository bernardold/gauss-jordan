#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "gauss-jordan.h"
#include "column.h"

extern gauss_jordan gj;

void init(int augmented_n, float** augmented_m, int groupDistribution) {
    gj.dimension = augmented_n - 1;
    gj.groupDistribution = groupDistribution;
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
    // Init MPI
    MPI_Init(NULL, NULL);
    // Initialize global Gauss Jordan info struct
    MPI_Comm_size(MPI_COMM_WORLD, &(gj.proc_num));
    gj.groupNumber = gj.dimension / gj.proc_num;
    
    // Send columns
    // Populate yours
    // Wait on sent data
}

int is_my_column(int col_idx, int my_rank) {
    return process_of_column(col_idx) == my_rank;
}

int process_of_column(int col_idx) {
    if (gj.groupDistribution) {
        return col_idx / gj.groupNumber;
    }
    else {
        return col_idx % gj.proc_num;
    }
}

