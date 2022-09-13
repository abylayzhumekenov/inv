#include <time.h>
#include <math.h>

#include <pcg_variants.h>
#include <petsc.h>
#include <petscblaslapack.h> 
#include "/home/abylay/KAUST/Libraries/C/petsc/src/ksp/ksp/impls/gmres/gmresimpl.h"

#include "inv_sampler.h"
#include "inv_random.h"

#define INV_SAMPLER_VERBOSE "SAMPLER:\t"


KSP sampler_ksp(MPI_Comm comm, Mat Q, int max_niter, int verbose){

    /* Create a GMRES solver */
    if(verbose) printf("%sSetting up a solver..\n", INV_SAMPLER_VERBOSE);
    KSP ksp;
    PC pc;
    KSPCreate(comm, &ksp);
    KSPSetOperators(ksp, Q, Q);
    KSPSetType(ksp, KSPGMRES);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCICC);  // Later fix to some parallel preconditioner
    PCSetFromOptions(pc);
    KSPSetPCSide(ksp, PC_SYMMETRIC);
    KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, max_niter);
    KSPSetComputeEigenvalues(ksp, PETSC_TRUE);
    KSPSetComputeRitz(ksp, PETSC_TRUE);
    KSPGMRESSetRestart(ksp, max_niter);
    KSPGMRESSetOrthogonalization(ksp,KSPGMRESClassicalGramSchmidtOrthogonalization);
    KSPGMRESSetCGSRefinementType(ksp, KSP_GMRES_CGS_REFINE_IFNEEDED);                   // Without this, even CGS fails...
    // KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization);  // MGS fails too...
    KSPSetFromOptions(ksp);
    if(verbose) printf("%sSolver created.\n", INV_SAMPLER_VERBOSE);

    return ksp;
}


Vec* sampler_gmrf(KSP ksp, Vec* z, int n_samples, int verbose, int verbose_sampler){

    /* Allocate space for x */
    if(verbose) printf("%sAllocating memory..\n", INV_SAMPLER_VERBOSE);
    Vec* x = (Vec*)malloc(n_samples * sizeof(Vec));    // Do it another way?

    for(int k=0; k<n_samples; k++){

        if(verbose && verbose_sampler) printf("%sSample #%i..\n", INV_SAMPLER_VERBOSE, k);

        /* Solve a system */
        if(verbose && verbose_sampler) printf("%sSolving a system..\n", INV_SAMPLER_VERBOSE);
        int niter, max_niter;
        Vec y;
        VecDuplicate(z[0], &y);
        VecDuplicate(z[k], &x[k]);
        KSPSolve(ksp, z[k], y);
        KSPGetIterationNumber(ksp, &niter);
        KSP_GMRES* gmres = (KSP_GMRES*)ksp->data;
        max_niter = gmres->max_k;
        if(verbose && verbose_sampler) printf("%sSolved in %i iterations.\n", INV_SAMPLER_VERBOSE, niter);

        /* Get Krylov vectors and the Hessenberg matrix */
        if(verbose && verbose_sampler) printf("%sCopying GMRES data..\n", INV_SAMPLER_VERBOSE);
        Vec* v = gmres->vecs + 2;
        double eig_diag[niter], eig_offdiag[niter-1];
        for(int i=0; i<niter; i++){
            eig_diag[i] = gmres->hes_origin[(max_niter+1)*i+i];
            if(i<niter-1) eig_offdiag[i] = gmres->hes_origin[(max_niter+1)*i+i+1];
        }

        /* Eigendecomposition of the Hessenberg matrix */
        if(verbose && verbose_sampler) printf("%sEigendecomposition..\n", INV_SAMPLER_VERBOSE);
        const char eig_v = 'I';
        int eig_n = niter, eig_info;
        double eig_work[5*niter];
        double* eig_vectors = (double*)malloc(niter*niter * sizeof(double));
        dsteqr_(&eig_v, &eig_n, eig_diag, eig_offdiag, eig_vectors, &eig_n, eig_work, &eig_info);
        if(verbose && verbose_sampler) printf("%sEigensolve finished.\n", INV_SAMPLER_VERBOSE);

        /* Compute right term of x = V U D^-1/2 U^T (beta e_1) = V * udube, that is, of 'udube'     }:8)     */
        if(verbose && verbose_sampler) printf("%sComputing weights..\n", INV_SAMPLER_VERBOSE);
        double beta0 = gmres->rnorm0;
        double udube[niter];
        for(int i=0; i<niter; i++){
            udube[i] = 0;
            for(int j=0; j<niter; j++){
                udube[i] += beta0 / sqrt(eig_diag[j]) * eig_vectors[j*niter] * eig_vectors[j*niter+i];
            }
        }

        /* Compute the sample x = V * udube */
        if(verbose && verbose_sampler) printf("%sAssembling a sample..\n", INV_SAMPLER_VERBOSE);
        VecSet(x[k], 0.0);
        VecAssemblyBegin(x[k]);
        VecAssemblyEnd(x[k]);
        VecMAXPY(x[k], niter, udube, v);

        /* Clean up */
        VecDestroy(&y);
        free(eig_vectors);
        if(verbose && verbose_sampler) printf("%sSample #%i finished.\n", INV_SAMPLER_VERBOSE, k);
    }

    if(verbose) printf("%sSuccessfully finished.\n", INV_SAMPLER_VERBOSE);

    return x;
}


Vec* sampler_std_normal(MPI_Comm comm, int n, int n_samples, int verbose){

    /* Create a pseudo-random number generator */
    if(verbose) printf("%sSetting seed..\n", INV_SAMPLER_VERBOSE);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    pcg64_random_t rng;
    pcg64_srandom_r(&rng, time(NULL), 0);   // need different seeds
    // pcg64_srandom_r(&rng, time(NULL), (intptr_t)rng);  // no need to set different seeds for different processes?

    // std::random_device sample_rd{};
    // std::mt19937 sample_gen{sample_rd()};
    // sample_gen.seed(time(0) * (rank + 1));    // use different seed for different processes, otherwise they produce the same numbers
    // std::normal_distribution<> sample_d{0,1};

    /* Generate random numbers in parallel */
    if(verbose) printf("%sAllocating memory..\n", INV_SAMPLER_VERBOSE);
    Vec* z = (Vec*)malloc(n_samples * sizeof(Vec));
    int Istart, Iend;
    VecCreate(comm, &z[0]);
    VecSetType(z[0], VECMPI);
    VecSetSizes(z[0], PETSC_DECIDE, n);
    VecSetFromOptions(z[0]);
    VecSetUp(z[0]);
    VecGetOwnershipRange(z[0], &Istart, &Iend);
    // VecDuplicateVecs(z[0], n_samples-1, &z);

    if(verbose) printf("%sGenerating standard normal..\n", INV_SAMPLER_VERBOSE);
    for(int i=0; i<n_samples; i++){
        if(i) VecDuplicate(z[0], &z[i]);
        for(int j=Istart; j<Iend; j++){
            VecSetValue(z[i], j, random_std_normal(&rng), INSERT_VALUES);
            // VecSetValue(z[i], j, sample_d(sample_gen), INSERT_VALUES);
        }
        VecAssemblyBegin(z[i]);
        VecAssemblyEnd(z[i]);
    }

    return z;
}