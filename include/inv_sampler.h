#ifndef INV_SAMPLER_H
#define INV_SAMPLER_H

#include <petsc.h>
#include "inv_random.h"


// KSP InvSamplerCreateKSP(MPI_Comm comm, Mat Q, int max_niter, int verbose);


PetscErrorCode InvSamplerStdNormal(pcg64_random_t* rng, Vec* z, int verbose);


PetscErrorCode InvSamplerGMRF(KSP ksp, Vec z, Vec* x, int verbose);


#endif