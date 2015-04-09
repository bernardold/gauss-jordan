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

    typedef struct {
        int groupDistribution;
        int proc_num;
        int groupNumber;
        int dimension;
    } gauss_jordan;

    void init(int augmented_n, float** augmented_m, int groupDistribution);
    int is_my_column(int col_idx, int my_rank);
    int process_of_column(int col_idx);
    

#ifdef	__cplusplus
}
#endif

#endif	/* GAUSS_JORDAN_H */

