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

Compile by running `make` in the project directory. Compilation is the same on PC, server or Ibex cluster.

## Dependencies

* petsc (configure with `--download-mumps --download-scalapack --download-metis --download-parmetis`)
* pcg_random
* lapack (shipped with petsc)
* mpi (shipped with petsc)
* mumps (shipped with petsc)
* metis (shipped with petsc)
* gklib (if you build metis yourself)

## Run with options

Run as `mpiexec [mpi_options] ./bin/main [inv_options]`

MPI options:
* `-n` number of processes

Program options:
* `-ns` number of samples
* `-nn` number of neighbors
* `-ni` max number of Krylov iterations
* `-tauy` precision for observations
* `-taub` precision for fixed effects
* `-v` verbose
* `-p` basic time and memory profiling

## Input and output

Run `Rscript generate.R [rscript_options]` from inside the `rscript/simulation` folder to generate precision matrices for simulated example, or from `rscript/application` folder to generate them for application to US temperature data. For the simulation, the hyperparameters can be set through command line options as:

* `-ms` space resolution
* `-mt` time resolution
* `-rs` space range
* `-rt` time range
* `-ss` marginal variance

This will generate the precision matrices `J.0`, `J.1`, `J.2`, `K.3`, `K.2`, `K.1`, projection matrices `A.t`, `A.s`, covariates `A.b` and observations `y`

The output will be written as binary vectors `d` (variance) and `mu` (mean), where the first 64 bits (PETSc header) can be ignored

Note that on linux machines PETSc uses a swapped endian, this must be taken into account when reading `d` and `mu`
