#include <gmresimpl.h>
#include <petscblaslapack.h> 

#include "inv_sampler.h"
#include "inv_shell.h"


PetscErrorCode InvSamplerStdNormal(pcg64_random_t* rng, Vec* z){

    int i_start, i_end;
    VecGetOwnershipRange(*z, &i_start, &i_end);
    for(int i=i_start; i<i_end; i++) VecSetValue(*z, i, InvRandomStdNormal(rng), INSERT_VALUES);
    VecAssemblyBegin(*z);
    VecAssemblyEnd(*z);

    return 0;
}


PetscErrorCode InvSamplerGMRF(KSP ksp, Vec z, Vec* x){
        double tt0, tt1, tt2, tt3, tt4, tt5, tt6, tt7;
        PetscTime(&tt0);
    /* Solve a system */
    int niter, max_niter;
    Vec xx, y;
    VecDuplicate(z, &y);
    VecDuplicate(z, &xx);
    VecSet(xx, 0.0);
    VecAssemblyBegin(xx);
    VecAssemblyEnd(xx);
    KSPSolve(ksp, z, y);
    KSPGetIterationNumber(ksp, &niter);
    KSP_GMRES* gmres = (KSP_GMRES*)ksp->data;
    max_niter = gmres->max_k;
        PetscTime(&tt1);

    /* Get Krylov vectors and the Hessenberg matrix */
    Vec* v = gmres->vecs + 2;
    double eig_diag[niter], eig_offdiag[niter-1];
    for(int i=0; i<niter; i++){
        eig_diag[i] = gmres->hes_origin[(max_niter+1)*i+i];
        if(i<niter-1) eig_offdiag[i] = gmres->hes_origin[(max_niter+1)*i+i+1];
    }
        PetscTime(&tt2);

    /* Eigendecomposition of the Hessenberg matrix */
    const char eig_v = 'I';
    int eig_n = niter, eig_info;
    double eig_work[5*niter];
    double* eig_vectors = (double*)malloc(niter*niter * sizeof(double));
    dsteqr_(&eig_v, &eig_n, eig_diag, eig_offdiag, eig_vectors, &eig_n, eig_work, &eig_info);
        PetscTime(&tt3);

    /* Compute right term of xx = V U D^-1/2 U^T (beta e_1) = V * udube, that is, of 'udube'     }:8)     */
    double beta0 = gmres->rnorm0;
    double udube[niter];
    for(int i=0; i<niter; i++){
        udube[i] = 0;
        for(int j=0; j<niter; j++){
            udube[i] += beta0 / sqrt(eig_diag[j]) * eig_vectors[j*niter] * eig_vectors[j*niter+i];
        }
    }
        PetscTime(&tt4);

    /* Compute the sample xx = V * udube */
    VecMAXPY(xx, niter, udube, v);
        PetscTime(&tt5);

    /* Unwind the preconditioner x = L^-T * V * xx */
    PC pc;
    KSPGetPC(ksp, &pc);
    PCApplySymmetricRight(pc, xx, *x);
        PetscTime(&tt6);

    /* Clean up */
    VecDestroy(&y);
    VecDestroy(&xx);
    free(eig_vectors);
        PetscTime(&tt7);

        tt7 = tt7 - tt6;
        tt6 = tt6 - tt5;
        tt5 = tt5 - tt4;
        tt4 = tt4 - tt3;
        tt3 = tt3 - tt2;
        tt2 = tt2 - tt1;
        tt1 = tt1 - tt0;
        tt0 = tt1 + tt2 + tt3 + tt4 + tt5 + tt6 + tt7;

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank) printf("\n\t\tSampling:\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", tt1/tt0, tt2/tt0, tt3/tt0, tt4/tt0, tt5/tt0, tt6/tt0, tt7/tt0);

    return 0;
}