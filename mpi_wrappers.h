/* 
 * File:   mpi_wrappers.h
 * Author: chris
 *
 * Created on April 9, 2015, 2:42 PM
 */

#ifndef MPI_WRAPPERS_H
#define	MPI_WRAPPERS_H

#ifdef	__cplusplus
extern "C" {
#endif

    void wait_all_wrapper(MPI_Request* requests, int size);
    void wait_wrapper(MPI_Request* request);

#ifdef	__cplusplus
}
#endif

#endif	/* MPI_WRAPPERS_H */

