#include <stdio.h>
#include "column.h"
#include "gauss-jordan.h"

extern gauss_jordan gj;

void print_column(float* data) {
    printf("Column: ");
    int i;
    for (i = 0; i < gj.dimension; i++) {
        printf("%f | ", data[i]);
    }
    printf("\n");
}