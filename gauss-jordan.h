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

    column** init(int augmented_n, float** augmented_m, int groupDistribution);
    

#ifdef	__cplusplus
}
#endif

#endif	/* GAUSS_JORDAN_H */

