#ifndef INV_SAMPLER_H
#define INV_SAMPLER_H

#include <petsc.h>


/**
 * @brief Create a KSP solver for sampling
 * 
 * @param comm 
 * @param Q 
 * @param max_niter 
 * @param verbose 
 * @return KSP 
 */
KSP sampler_ksp(MPI_Comm comm, Mat Q, int max_niter, int verbose);


/**
 * @brief Sample from a GMRF by solving Qx = z, using KSP
 * 
 * @param ksp 
 * @param z 
 * @param n_samples 
 * @param verbose 
 * @return Vec* 
 */
Vec* sampler_gmrf(KSP ksp, Vec* z, int n_samples, int verbose, int verbose_sampler);


/**
 * @brief Sample from a standard normal distribution
 * 
 * @param comm 
 * @param n 
 * @param n_samples 
 * @param verbose 
 * @return Vec* 
 */
Vec* sampler_std_normal(MPI_Comm comm, int n, int n_samples, int verbose);

#endif