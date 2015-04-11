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

    struct Column_Type {
        float* data;
        int idx;
    };
    typedef struct Column_Type* column_t;
    
    column_t create_column(int idx, float* data);
    void delete_column(column_t col_ptr);
    void print_column(column_t column);


#ifdef	__cplusplus
}
#endif

#endif	/* COLUMN_H */

