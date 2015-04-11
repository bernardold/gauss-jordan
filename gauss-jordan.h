/* 
 * Project: Gauss-Jordan elimination method by columns (KJI-form)
 *          Implemented with MPI
 * File:   gauss-jordan.h
 * Author: Chris Aslanoglou
 */

#ifndef GAUSS_JORDAN_H
#define	GAUSS_JORDAN_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "column.h"

    typedef struct {
        int groupDistribution;
        int proc_num;
        int groupNumber;
        int dimension;
    } gauss_jordan;

    column_t* init(int augmented_n, float** augmented_m, int groupDistribution);
    void destroy(column_t** my_cols_ptr, int my_rank);
    

#ifdef	__cplusplus
}
#endif

#endif	/* GAUSS_JORDAN_H */

