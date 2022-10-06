# inv

C library for sparse matrix inversion

Currently computes only the diagonal of the inverse

## Compile

Compile by running `make` in the project directory

## Dependencies

* petsc (configure with `--download-mumps --download-scalapack`)
* mpi (loaded with petsc)
* metis
* lapack
* pcg_random

## Run with options

Run as `mpiexec [mpi_options] ./bin/main [inv_options]`

MPI options:
* `-n` number of processes

Program options:
* `-fin` MatrixMarket file name
* `-ns` number of samples
* `-nn` number of neighbors
* `-nmax` max number of iterations
* `-v` verbose
* `-vs` separate verbose for sampling phase

## Input and output

Place input matrix files `filename` in the `data` folder

The output will be written as a binary vector `x_[filename]`, where the first 64 bits (header) can be ignored
