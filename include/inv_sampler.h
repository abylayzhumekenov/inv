#ifndef INV_SAMPLER_H
#define INV_SAMPLER_H

#include <petsc.h>

#include "inv_random.h"


/**
 * @brief Sample from IID standard normal
 * 
 * @param rng 
 * @param z 
 * @return PetscErrorCode 
 */
PetscErrorCode InvSamplerStdNormal(pcg64_random_t* rng, Vec* z);


/**
 * @brief Sample from a GMRF using a Krylov solver and a standard normal vector
 * 
 * @param ksp 
 * @param z 
 * @param x 
 * @return PetscErrorCode 
 */
PetscErrorCode InvSamplerGMRF(KSP ksp, Vec z, Vec* x);


#endif