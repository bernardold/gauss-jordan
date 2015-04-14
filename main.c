/* 
 * Project: Gauss-Jordan elimination method by columns (KJI-form)
 *          Implemented with MPI
 * File:   main.c
 * Author: Chris Aslanoglou
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>

#include "gauss_jordan.h"

gauss_jordan gj;

int main(int argc, char** argv) {
    srand(time(NULL));
    if (argc != 3) {
        fprintf(stderr, "Usage: %s [dimension (not augmented)] [use k-group distribution]\n",
                argv[0]);
        exit(EXIT_SUCCESS);
    }
    int dimension = atoi(argv[1]);
    if (dimension <= 0 || (dimension % 2 != 0)) {
        fprintf(stderr, "[dimension] must be a positive even number\n");
        exit(EXIT_SUCCESS);
    }
    int groupsDistribution = atoi(argv[2]);
    if (groupsDistribution < 0 || groupsDistribution > 1) {
        fprintf(stderr, "[use k-group distribution must be either 0 or 1]\n");
        exit(EXIT_SUCCESS);
    }
    // Create augmented random array
    float** augmented_m = create_augmented_matrix(dimension);
    int augmented_n = dimension + 1;
    /* Initial test case
    int augmented_n = 5;
    float **augmented_m = NULL;
    augmented_m = malloc((augmented_n - 1) * sizeof(float*));
    int i;
    for (i = 0; i < augmented_n - 1; i++) {
        augmented_m[i] =  malloc(augmented_n * sizeof(float));
    }
    augmented_m[0][0] = 0;
    augmented_m[0][1] = 2;
    augmented_m[0][2] = 0;
    augmented_m[0][3] = 1;
    augmented_m[0][4] = 0;
    augmented_m[1][0] = 2;
    augmented_m[1][1] = 2;
    augmented_m[1][2] = 3;
    augmented_m[1][3] = 2;
    augmented_m[1][4] = -2;
    augmented_m[2][0] = 4;
    augmented_m[2][1] = -3;
    augmented_m[2][2] = 0;
    augmented_m[2][3] = 1;
    augmented_m[2][4] = -7;
    augmented_m[3][0] = 6;
    augmented_m[3][1] = 1;
    augmented_m[3][2] = -6;
    augmented_m[3][3] = -5;
    augmented_m[3][4] = 6;
     * */
    
    column_t* my_cols = init(augmented_n, augmented_m, groupsDistribution);
    MPI_Barrier(MPI_COMM_WORLD);
    gj_kgi_main_loop(my_cols);
    if (gj.my_rank == 0) {
        // Master process has b-vector
//        printf("Master printing solution\n");
//        print_column(my_cols[gj.group_number]);
    }
    // All processes must reach to the barrier for the de-allocation process and
    // the termination of the MPI framework, or else, master process doesn't get
    // finalize properly
    MPI_Barrier(MPI_COMM_WORLD);
    destroy(&my_cols);
    return (EXIT_SUCCESS);
}

