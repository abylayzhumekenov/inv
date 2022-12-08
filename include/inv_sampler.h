#ifndef INV_SAMPLER_H
#define INV_SAMPLER_H

#include <petsc.h>
#include "inv_random.h"
#include "inv_is.h"


KSP InvSamplerCreateKSP(MPI_Comm comm, Mat Q, int max_niter, int verbose);


Vec InvSamplerStdNormal(MPI_Comm comm, Vec z, pcg64_random_t* rng, InvIS* mapping, int verbose);


Vec InvSamplerGMRF(KSP ksp, Vec x, Vec z, int verbose);


#endif