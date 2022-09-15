#include <time.h>
#include <math.h>

#include <pcg_variants.h>
#include <petsc.h>
#include <petscblaslapack.h> 
#include <gmresimpl.h>
#include <sbaij.h>

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
    PCSetUp(pc);
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


Vec sampler_std_normal(MPI_Comm comm, pcg64_random_t* rng, int n, int verbose, int k){

    /* Generate random numbers in parallel */
    if(verbose) printf("%sSample #%i..\n", INV_SAMPLER_VERBOSE, k);
    if(verbose) printf("%sSampling standard normal..\n", INV_SAMPLER_VERBOSE);
    Vec z;
    int Istart, Iend;
    VecCreate(comm, &z);
    VecSetType(z, VECMPI);
    VecSetSizes(z, PETSC_DECIDE, n);
    VecSetFromOptions(z);
    VecSetUp(z);
    VecGetOwnershipRange(z, &Istart, &Iend);
    for(int j=Istart; j<Iend; j++){
        VecSetValue(z, j, random_std_normal(rng), INSERT_VALUES);
    }
    VecAssemblyBegin(z);
    VecAssemblyEnd(z);

    return z;
}


Vec sampler_premultiply_z(KSP ksp, Vec z, int verbose){

    /* 
    ICC is stored as M = D^-2 + U
    where D = diag(L)
    and U = D^-1 (-L^T + D)
    thus, L = (I - U^T) D

    Note: D is at the end of array

    ->mbs = int (or n) (number of rows/columns)
    ->nz = int (or nz) (number of nonzeros)
    ->a = double* of size nz (stores values)
    ->i = int* of size n+1 (stores ranges of each row)
    ->j = int* of size nz (stores column indices, diagonals at the end)
    ->diag = int* of size n (stores indices to diagonals in array j, e.g. j[diag])      */

    if(verbose) printf("%sPrepmultiplying standard normal..\n", INV_SAMPLER_VERBOSE);
    Mat L;
    PC pc;
    KSPGetPC(ksp, &pc);
    PCFactorGetMatrix(pc, &L);
    Mat_SeqSBAIJ *LL = (Mat_SeqSBAIJ*)L->data;
    double *z_array, zz_array[LL->mbs];
    VecGetArray(z, &z_array);

    for(int i=0; i<LL->mbs; i++){
        zz_array[i] = z_array[i] / sqrt(LL->a[LL->diag[i]]);
        z_array[i] = zz_array[i];
    }
    for(int i=0; i<LL->mbs; i++){
        for(int j=LL->i[i]; j<LL->i[i+1]-1; j++){
            z_array[LL->j[j]] -= LL->a[j] * zz_array[i];
        }
    }
    VecRestoreArray(z, &z_array);
    if(verbose) printf("%sPrepmultiplying done.\n", INV_SAMPLER_VERBOSE);

    return z;
}


Vec sampler_gmrf(KSP ksp, Vec z, int verbose, int k){

    if(verbose) printf("%sSampling GMRF..\n", INV_SAMPLER_VERBOSE);

    /* Solve a system */
    if(verbose) printf("%sSolving a system..\n", INV_SAMPLER_VERBOSE);
    int niter, max_niter;
    Vec x, xx, y;
    VecDuplicate(z, &y);
    VecDuplicate(z, &x);
    VecDuplicate(z, &xx);
    VecSet(x, 0.0);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    z = sampler_premultiply_z(ksp, z, verbose);
    KSPSolve(ksp, z, y);
    KSPGetIterationNumber(ksp, &niter);
    KSP_GMRES* gmres = (KSP_GMRES*)ksp->data;
    max_niter = gmres->max_k;
    if(verbose) printf("%sSolved in %i iterations.\n", INV_SAMPLER_VERBOSE, niter);

    /* Get Krylov vectors and the Hessenberg matrix */
    if(verbose) printf("%sCopying GMRES data..\n", INV_SAMPLER_VERBOSE);
    Vec* v = gmres->vecs + 2;
    double eig_diag[niter], eig_offdiag[niter-1];
    for(int i=0; i<niter; i++){
        eig_diag[i] = gmres->hes_origin[(max_niter+1)*i+i];
        if(i<niter-1) eig_offdiag[i] = gmres->hes_origin[(max_niter+1)*i+i+1];
    }

    /* Eigendecomposition of the Hessenberg matrix */
    if(verbose) printf("%sEigendecomposition..\n", INV_SAMPLER_VERBOSE);
    const char eig_v = 'I';
    int eig_n = niter, eig_info;
    double eig_work[5*niter];
    double* eig_vectors = (double*)malloc(niter*niter * sizeof(double));
    dsteqr_(&eig_v, &eig_n, eig_diag, eig_offdiag, eig_vectors, &eig_n, eig_work, &eig_info);
    if(verbose) printf("%sEigensolve finished.\n", INV_SAMPLER_VERBOSE);

    /* Compute right term of x = V U D^-1/2 U^T (beta e_1) = V * udube, that is, of 'udube'     }:8)     */
    if(verbose) printf("%sComputing weights..\n", INV_SAMPLER_VERBOSE);
    double beta0 = gmres->rnorm0;
    double udube[niter];
    for(int i=0; i<niter; i++){
        udube[i] = 0;
        for(int j=0; j<niter; j++){
            udube[i] += beta0 / sqrt(eig_diag[j]) * eig_vectors[j*niter] * eig_vectors[j*niter+i];
        }
    }

    /* Compute the sample x = V * udube */
    if(verbose) printf("%sAssembling a sample..\n", INV_SAMPLER_VERBOSE);
    VecMAXPY(x, niter, udube, v);

    /* Unwind the preconditioner x = L^-T * V */
    if(verbose) printf("%sUnwinding preconditioner..\n", INV_SAMPLER_VERBOSE);
    PC pc;
    KSPGetPC(ksp, &pc);
    PCApplySymmetricRight(pc, x, xx);

    /* Clean up */
    VecDestroy(&y);
    VecDestroy(&x);
    free(eig_vectors);
    if(verbose) printf("%sSample #%i finished.\n", INV_SAMPLER_VERBOSE, k);

    return xx;
}
