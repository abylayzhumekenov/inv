# inv

C program for selected inversion of sparse precision matrices in latent Gaussian models

The latent Gaussian model are hierarchical Bayesian models

$$
    something
$$


## Compile

Compile by running `make` in the project directory

Run `make server` to compile on a server and `make cluster` to compile on Ibex cluster

## Dependencies

* petsc (configure with `--download-mumps --download-scalapack`)
* mpi (loaded with petsc)
* metis (can use the one shipped with petsc)
* gklib (if you build metis yourself)
* lapack
* pcg_random

## Run with options

Run as `mpiexec [mpi_options] ./bin/main [inv_options]`

MPI options:
* `-n` number of processes

Program options:
* `-ns` number of samples
* `-nn` number of neighbors
* `-nmax` max number of iterations
* `-sel` selected inversion
* `-tauy` precision for observations
* `-taub` precision for fixed effects
* `-v` verbose
* `-vs` separate verbose for sampling phase

## Input and output

Run `Rscript generate_Q.R` from inside the `data` folder to generate precision matrices with parameters

* `-ms` space resolution
* `-mt` time resolution
* `-rt` time range

This will generate the precision matrices `J.0`, `J.1`, `J.2`, `K.3`, `K.2`, `K.1` and observations `y`

The output will be written as binary vectors `d` (variance) and `mu` (mean), where the first 64 bits (PETSc header) can be ignored

Note that on linux machines PETSc uses a swapped endian, this must be taken into account when reading `d` and `mu`
