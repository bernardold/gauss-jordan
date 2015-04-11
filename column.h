/* 
 * Project: Gauss-Jordan elimination method by columns (KJI-form)
 *          Implemented with MPI
 * File:   column.h
 * Author: Chris Aslanoglou
 */

#ifndef COLUMN_H
#define	COLUMN_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct Column_Type {
        float* data;
        int idx;
    } column;
    
    column* create_column(int idx, float* data);
    void delete_column(column** col_ptr);
    void print_column(column* column);


#ifdef	__cplusplus
}
#endif

#endif	/* COLUMN_H */

