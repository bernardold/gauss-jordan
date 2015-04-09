/* 
 * Project: Gauss-Jordan elimination method by columns (KJI-form)
 *          Implemented with MPI
 * File:   main.c
 * Author: Chris Aslanoglou
 */

#include <stdio.h>
#include <stdlib.h>

#include "gauss-jordan.h"

gauss_jordan gj;

int main(int argc, char** argv) {
    // WIP - Fixed params
    int groupsDistribution = 1;
    // Read array
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
    
    init(augmented_n, augmented_m, groupsDistribution);
    

    return (EXIT_SUCCESS);
}
