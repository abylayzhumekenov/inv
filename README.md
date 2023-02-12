# inv

C program for selected inversion of sparse precision matrices in latent Gaussian models.

The latent Gaussian model are hierarchical Bayesian models. In our case, we solve the Gaussian problem with fixed hyperparameters:

$$
\begin{aligned}
    y|x,\theta &\sim \mathcal{N}(Ax,Q_y^{-1}) \\
    x|\theta &\sim \mathcal{N}(0,Q_x^{-1}) \\
    \theta &= \theta_0 (\text{fixed})
\end{aligned}
$$

The posterior mean and marginal variances are

$$
    \mu = Q^{-1}b \quad\text{and}\quad \text{diag}\Sigma = \text{diag}Q^{-1}
$$

where $Q=Q_x+A^TQ_yA$ is a sparse posterior precision matrix.


## Compile

Compile by running `make` in the project directory, run `make server` to compile on a server or Ibex cluster.

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
* `-prof` basic time and memory profiling

## Input and output

Run `Rscript generate_Q.R` from inside the `data` folder to generate precision matrices with parameters

* `-ms` space resolution
* `-mt` time resolution
* `-rt` time range

This will generate the precision matrices `J.0`, `J.1`, `J.2`, `K.3`, `K.2`, `K.1` and observations `y`

The output will be written as binary vectors `d` (variance) and `mu` (mean), where the first 64 bits (PETSc header) can be ignored

Note that on linux machines PETSc uses a swapped endian, this must be taken into account when reading `d` and `mu`
