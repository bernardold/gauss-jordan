Gauss Jordan Elimination - KJI form
=======

The program implements the Gauss Jordan Elimination method for solving a linear
system of equations, more on the method [here](http://en.wikipedia.org/wiki/Gaussian_elimination).

The program is implemented in C, using MPI [v. 3.0.4] and uses the KJI form of the algorithm, which means elimination by columns.

##Notes on the parallelism of the problem
- The augmented matrix (A|b) is needed for the initiallization
- The columns of the matrix are distributed in two ways, either in k-groups (that is k-consecutive columns
are given to each process in a sequential fashion), or with circular-distribution (card-shuffling fashion).
- The b-vector is handled by the master process at all times

##Misc notes
- In the initial version, the NxN array is generated randomly (*future work:* reading the array from a csv file)
- A makefile and a test file are included
- The project was developped for a university project for the Parallel Algorithms course.
