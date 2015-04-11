#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "column.h"
#include "gauss-jordan.h"

extern gauss_jordan gj;

/**
 * Creates a Column datatype (allocates the memory
 * needed and returns a pointer to it)
 * @param idx The index of the column in the original matrix
 * @param data The coefficients of the column
 * @return A pointer to the allocated Column datatype
 */
column_t create_column(int idx, float* data) {
    column_t col = malloc(sizeof(struct Column_Type));
    col->idx = idx;
    col->data = malloc(gj.dimension * sizeof(float));
    memcpy(col->data, data, gj.dimension * sizeof(float));
    return col;
}

/**
 * Deletes the space allocated for the column and nullifies the
 * pointer to it
 * @param col_ptr The address of the pointer to the column datatype
 */
void delete_column(column_t col_ptr) {
    (col_ptr)->idx = -1;
    free((col_ptr)->data);
    (col_ptr)->data = NULL;
}

/**
 * Prints the column of length dimension (defined by the global 
 * Gauss Jordan struct)
 * @param column A pointer to the column to be printed
 */
void print_column(column_t column) {
    if (column->data != NULL) {
        printf("Column[%d]: ", column->idx);
        int i;
        for (i = 0; i < gj.dimension; i++) {
            printf("%f | ", column->data[i]);
        }
        printf("\n");
    }
    else {
        printf("Column isn't instantiated\n");
    }
}