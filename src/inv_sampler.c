#include <gmresimpl.h>
#include <petscblaslapack.h> 

#include "inv_sampler.h"
#include "inv_shell.h"

#define INV_SAMPLER_VERBOSE "SAMPLER:\t"


KSP InvSamplerCreateKSP(MPI_Comm comm, Mat Q, int max_niter, int verbose){

    if(verbose) printf("%sSetting up a solver with a shell preconditioner...\n", INV_SAMPLER_VERBOSE);
    InvShellPC* shell;
    KSP ksp;
    PC pc;
    
    KSPCreate(comm, &ksp);
    KSPSetOperators(ksp, Q, Q);
    KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, max_niter);
    KSPSetComputeEigenvalues(ksp, PETSC_TRUE);
    KSPSetPCSide(ksp, PC_SYMMETRIC);
    KSPGetPC(ksp, &pc);

    PCSetType(pc, PCSHELL);
    InvShellCreate(&shell);
    PCShellSetApply(pc, InvShellApply);
    PCShellSetApplyTranspose(pc, InvShellApplyTranspose);
    PCShellSetApplyBA(pc, InvShellApplyBA);
    PCShellSetApplySymmetricLeft(pc, InvShellApplyLeft);
    PCShellSetApplySymmetricRight(pc, InvShellApplyRight);
    PCShellSetContext(pc, shell);
    PCShellSetDestroy(pc, InvShellDestroy);
    PCShellSetName(pc, "shell");
    InvShellSetup(pc, Q, comm);
    if(verbose) printf("%sSolver created.\n", INV_SAMPLER_VERBOSE);

    return ksp;
}


Vec InvSamplerStdNormal(MPI_Comm comm, Vec z, pcg64_random_t* rng, InvIS* mapping, int verbose){

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if(verbose) printf("%sSampling standard normal..\n", INV_SAMPLER_VERBOSE);
    int Istart, Iend;
    if(!z){
        VecCreate(comm, &z);
        VecSetType(z, VECMPI);
        VecSetSizes(z, mapping->n_local[rank], PETSC_DETERMINE);
        VecSetFromOptions(z);
        VecSetUp(z);
    }
    VecGetOwnershipRange(z, &Istart, &Iend);
    for(int j=Istart; j<Iend; j++){
        VecSetValue(z, j, InvRandomStdNormal(rng), INSERT_VALUES);
    }
    VecAssemblyBegin(z);
    VecAssemblyEnd(z);

    return z;
}


Vec InvSamplerGMRF(KSP ksp, Vec xx, Vec z, int verbose){

    if(verbose) printf("%sSampling GMRF..\n", INV_SAMPLER_VERBOSE);

    /* Solve a system */
    if(verbose) printf("%sSolving a system...\n", INV_SAMPLER_VERBOSE);
    int niter, max_niter;
    Vec x, y;
    VecDuplicate(z, &y);
    VecDuplicate(z, &x);
    if(!xx) VecDuplicate(z, &xx);
    VecSet(x, 0.0);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    KSPSolve(ksp, z, y);
    KSPGetIterationNumber(ksp, &niter);
    KSP_GMRES* gmres = (KSP_GMRES*)ksp->data;
    max_niter = gmres->max_k;
    if(verbose) printf("%sSolved in %i iterations.\n", INV_SAMPLER_VERBOSE, niter);

    /* Get Krylov vectors and the Hessenberg matrix */
    if(verbose) printf("%sCopying GMRES data...\n", INV_SAMPLER_VERBOSE);
    Vec* v = gmres->vecs + 2;
    double eig_diag[niter], eig_offdiag[niter-1];
    for(int i=0; i<niter; i++){
        eig_diag[i] = gmres->hes_origin[(max_niter+1)*i+i];
        if(i<niter-1) eig_offdiag[i] = gmres->hes_origin[(max_niter+1)*i+i+1];
    }

    /* Eigendecomposition of the Hessenberg matrix */
    if(verbose) printf("%sEigendecomposition...\n", INV_SAMPLER_VERBOSE);
    const char eig_v = 'I';
    int eig_n = niter, eig_info;
    double eig_work[5*niter];
    double* eig_vectors = (double*)malloc(niter*niter * sizeof(double));
    dsteqr_(&eig_v, &eig_n, eig_diag, eig_offdiag, eig_vectors, &eig_n, eig_work, &eig_info);
    if(verbose) printf("%sEigensolve finished.\n", INV_SAMPLER_VERBOSE);

    /* Compute right term of x = V U D^-1/2 U^T (beta e_1) = V * udube, that is, of 'udube'     }:8)     */
    if(verbose) printf("%sComputing weights...\n", INV_SAMPLER_VERBOSE);
    double beta0 = gmres->rnorm0;
    double udube[niter];
    for(int i=0; i<niter; i++){
        udube[i] = 0;
        for(int j=0; j<niter; j++){
            udube[i] += beta0 / sqrt(eig_diag[j]) * eig_vectors[j*niter] * eig_vectors[j*niter+i];
        }
    }

    /* Compute the sample x = V * udube */
    if(verbose) printf("%sAssembling a sample...\n", INV_SAMPLER_VERBOSE);
    VecMAXPY(x, niter, udube, v);

    /* Unwind the preconditioner x = L^-T * V */
    if(verbose) printf("%sUnwinding preconditioner...\n", INV_SAMPLER_VERBOSE);
    PC pc;
    KSPGetPC(ksp, &pc);
    PCApplySymmetricRight(pc, x, xx);

    /* Clean up */
    VecDestroy(&y);
    VecDestroy(&x);
    // free(eig_vectors);
    if(verbose) printf("%sSample computed.\n", INV_SAMPLER_VERBOSE);

    return xx;
}