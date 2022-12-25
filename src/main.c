#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <petsc.h>
#include <mpi.h>

#include "inv_matrix.h"
#include "inv_graph.h"
#include "inv_shell.h"
#include "inv_sampler.h"
#include "inv_random.h"


int main(int argc, char **argv){


    /* Set parameters from the options */
    int max_niter = 1000, n_sample = 1000, n_neighbor = 1, verbose = 0, verbose_s = 0, n_part = 1;
    for(int i=0; i<argc; i++){
        if(!strcmp(argv[i], "-ns")){
            n_sample = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nmax")){
            max_niter = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nn")){
            n_neighbor = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-v")){
            verbose = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-vs")){
            verbose_s = atoi(argv[i+1]);
        }
        argv[i] = 0;
    }


    /* Initialize MPI and PETSc */
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    PetscInitialize(&argc, &argv, NULL, NULL);
    verbose = (!rank) && verbose;


    /* Binary matrix load */
    PetscViewer viewer;
    Mat J0_global, J1_global, J2_global, J0, J1, J2, K1, K2, K3;
    MatCreate(PETSC_COMM_WORLD, &J0_global);
    MatCreate(PETSC_COMM_WORLD, &J1_global);
    MatCreate(PETSC_COMM_WORLD, &J2_global);
    MatCreate(PETSC_COMM_SELF, &J0);
    MatCreate(PETSC_COMM_SELF, &J1);
    MatCreate(PETSC_COMM_SELF, &J2);
    MatCreate(PETSC_COMM_SELF, &K1);
    MatCreate(PETSC_COMM_SELF, &K2);
    MatCreate(PETSC_COMM_SELF, &K3);
    MatSetType(J0, MATSEQAIJ);
    MatSetType(J1, MATSEQAIJ);
    MatSetType(J2, MATSEQAIJ);
    MatSetType(K1, MATSEQAIJ);
    MatSetType(K2, MATSEQAIJ);
    MatSetType(K3, MATSEQAIJ);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/J0", FILE_MODE_READ, &viewer);
    MatLoad(J0_global, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/J1", FILE_MODE_READ, &viewer);
    MatLoad(J1_global, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/J2", FILE_MODE_READ, &viewer);
    MatLoad(J2_global, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, "data/J0", FILE_MODE_READ, &viewer);
    MatLoad(J0, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, "data/J1", FILE_MODE_READ, &viewer);
    MatLoad(J1, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, "data/J2", FILE_MODE_READ, &viewer);
    MatLoad(J2, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, "data/K1", FILE_MODE_READ, &viewer);
    MatLoad(K1, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, "data/K2", FILE_MODE_READ, &viewer);
    MatLoad(K2, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, "data/K3", FILE_MODE_READ, &viewer);
    MatLoad(K3, viewer);


    /* Load and partition a temporal graph */
    InvGraph* graph = InvGraphCreate("data/J2", verbose);
    InvGraphPartition(graph, size, verbose);
    InvGraphNeighbor(graph, n_neighbor, verbose);
    // InvGraphPrint(graph, verbose);


    /* Get dimensions */
    int ns, nt_local, nt_ghost, nt_exten, nt_separ, nt, n_local, n_ghost, n_exten, n_separ, n;
    MatGetSize(K3, &ns, &ns);
    nt_local = graph->offset[1] - graph->offset[0];
    nt_ghost = graph->offset[2] - graph->offset[1];
    nt_exten = graph->offset[2] - graph->offset[0];
    nt_separ = graph->offset[3] - graph->offset[2];
    nt = graph->n_vert;
    n_local = nt_local * ns;
    n_ghost = nt_ghost * ns;
    n_exten = nt_exten * ns;
    n_separ = nt_separ * ns;
    n = ns * nt;


    /* Create index sets */
    IS ist_local, ist_ghost, ist_exten, ist_separ, is_local, is_ghost, is_exten, is_separ, is_origi;
    ISCreateGeneral(PETSC_COMM_SELF, nt_local, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &ist_local);
    ISCreateGeneral(PETSC_COMM_SELF, nt_ghost, graph->vert + graph->offset[1], PETSC_COPY_VALUES, &ist_ghost);
    ISCreateGeneral(PETSC_COMM_SELF, nt_exten, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &ist_exten);
    ISCreateGeneral(PETSC_COMM_SELF, nt_separ, graph->vert + graph->offset[2], PETSC_COPY_VALUES, &ist_separ);
    ISCreateBlock(PETSC_COMM_SELF, ns, nt_local, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &is_local);
    ISCreateBlock(PETSC_COMM_SELF, ns, nt_ghost, graph->vert + graph->offset[1], PETSC_COPY_VALUES, &is_ghost);
    ISCreateBlock(PETSC_COMM_SELF, ns, nt_exten, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &is_exten);
    ISCreateBlock(PETSC_COMM_SELF, ns, nt_separ, graph->vert + graph->offset[2], PETSC_COPY_VALUES, &is_separ);
    ISCreateStride(PETSC_COMM_SELF, n_local, 0, 1, &is_origi);

    
    // /* Create the local matrix */
    // Mat J0_local, J1_local, J2_local, Q_local, Q1_local, Q2_local, Q;
    // MatCreateSubMatrix(J0, ist_local, NULL, MAT_INITIAL_MATRIX, &J0_local);
    // MatCreateSubMatrix(J1, ist_local, NULL, MAT_INITIAL_MATRIX, &J1_local);
    // MatCreateSubMatrix(J2, ist_local, NULL, MAT_INITIAL_MATRIX, &J2_local);
    // MatSeqAIJKron(J0_local, K3, MAT_INITIAL_MATRIX, &Q_local);
    // MatSeqAIJKron(J1_local, K2, MAT_INITIAL_MATRIX, &Q1_local);
    // MatSeqAIJKron(J2_local, K1, MAT_INITIAL_MATRIX, &Q2_local);
    // MatAXPY(Q_local, 1.0, Q1_local, DIFFERENT_NONZERO_PATTERN);
    // MatAXPY(Q_local, 1.0, Q2_local, DIFFERENT_NONZERO_PATTERN);
    // MatCreateMPIMatConcatenateSeqMat(PETSC_COMM_WORLD, Q_local, n_local, MAT_INITIAL_MATRIX, &Q);

    Mat Q, Q1, Q2;
    MatMPIAIJKron(J0_global, K3, MAT_INITIAL_MATRIX, &Q);
    MatMPIAIJKron(J1_global, K2, MAT_INITIAL_MATRIX, &Q1);
    MatMPIAIJKron(J2_global, K1, MAT_INITIAL_MATRIX, &Q2);
    MatAXPY(Q, 1.0, Q1, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Q, 1.0, Q2, DIFFERENT_NONZERO_PATTERN);
    MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY);

    
    /* Create the extended matrix */
    Mat J0_exten, J1_exten, J2_exten, Q_exten, Q1_exten, Q2_exten;
    MatCreateSubMatrix(J0, ist_exten, ist_exten, MAT_INITIAL_MATRIX, &J0_exten);
    MatCreateSubMatrix(J1, ist_exten, ist_exten, MAT_INITIAL_MATRIX, &J1_exten);
    MatCreateSubMatrix(J2, ist_exten, ist_exten, MAT_INITIAL_MATRIX, &J2_exten);
    MatSeqAIJKron(J0_exten, K3, MAT_INITIAL_MATRIX, &Q_exten);
    MatSeqAIJKron(J1_exten, K2, MAT_INITIAL_MATRIX, &Q1_exten);
    MatSeqAIJKron(J2_exten, K1, MAT_INITIAL_MATRIX, &Q2_exten);
    MatAXPY(Q_exten, 1.0, Q1_exten, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Q_exten, 1.0, Q2_exten, DIFFERENT_NONZERO_PATTERN);
    MatAssemblyBegin(Q_exten, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q_exten, MAT_FINAL_ASSEMBLY);


    /* Create the separator portion of the matrix */
    Mat J0_separ, J1_separ, J2_separ, Q_separ, Q1_separ, Q2_separ;
    MatCreateSubMatrix(J0, ist_exten, ist_separ, MAT_INITIAL_MATRIX, &J0_separ);
    MatCreateSubMatrix(J1, ist_exten, ist_separ, MAT_INITIAL_MATRIX, &J1_separ);
    MatCreateSubMatrix(J2, ist_exten, ist_separ, MAT_INITIAL_MATRIX, &J2_separ);
    MatSeqAIJKron(J0_separ, K3, MAT_INITIAL_MATRIX, &Q_separ);
    MatSeqAIJKron(J1_separ, K2, MAT_INITIAL_MATRIX, &Q1_separ);
    MatSeqAIJKron(J2_separ, K1, MAT_INITIAL_MATRIX, &Q2_separ);
    MatAXPY(Q_separ, 1.0, Q1_separ, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Q_separ, 1.0, Q2_separ, DIFFERENT_NONZERO_PATTERN);
    MatAssemblyBegin(Q_separ, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q_separ, MAT_FINAL_ASSEMBLY);


    /* Compute the direct solution */   /* Ideally, should use Takahashi equations, but will not be able to solve with it... */
    Vec d;
    Mat L, II, S;
    IS is_direct;
    VecCreateSeq(PETSC_COMM_SELF, n_exten, &d);
    MatSetOption(Q_exten, MAT_SYMMETRIC, PETSC_TRUE);
    MatGetOrdering(Q_exten, MATORDERINGNATURAL, &is_direct, &is_direct);
    MatGetFactor(Q_exten, MATSOLVERMUMPS, MAT_FACTOR_CHOLESKY, &L);
    MatCholeskyFactorSymbolic(L, Q_exten, is_direct, NULL);
    MatCholeskyFactorNumeric(L, Q_exten, NULL);
    MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, n_exten, n_exten, NULL, &II);
    MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, n_exten, n_exten, NULL, &S);
    for(int i=0; i<n_exten; i++) MatSetValue(II, i, i, 1.0, INSERT_VALUES);
    MatAssemblyBegin(II, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(II, MAT_FINAL_ASSEMBLY);
    MatMatSolve(L, II, S);
    MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
    MatGetDiagonal(S, d);


    /* Create a solver for correction */
    KSP ksp_exten;
    PC pc_exten;
    KSPCreate(PETSC_COMM_SELF, &ksp_exten);
    KSPSetOperators(ksp_exten, Q_exten, Q_exten);
    KSPSetType(ksp_exten, KSPGMRES);
    KSPGetPC(ksp_exten, &pc_exten);
    PCSetType(pc_exten, PCBJACOBI);
    PCSetFromOptions(pc_exten);
    KSPSetTolerances(ksp_exten, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, max_niter);
    KSPSetFromOptions(ksp_exten);


    /* Create a solver for sampling */
    InvShellPC* shell;
    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
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
    InvShellSetup(pc, Q, PETSC_COMM_WORLD);


    /* Prepare sampling vectors */
    pcg64_random_t rng = InvRandomCreate(0, 0, 1);
    VecScatter scatter;
    Vec x, x_scatter, x_separ, z, y, w;
    int n_local_global;
    MatGetLocalSize(Q, &n_local_global, &n_local_global);
    VecCreateMPI(PETSC_COMM_WORLD, n_local_global, PETSC_DETERMINE, &x);
    VecDuplicate(x, &z);
    VecScatterCreateToAll(x, &scatter, &x_scatter);
    VecCreateSeq(PETSC_COMM_SELF, n_exten, &y);
    VecDuplicate(y, &w);


    /* Sample correction */
    for(int i=0; i<n_sample; i++){

        /* Sampling */
        InvSamplerStdNormal(&rng, &z, verbose_s);
        InvSamplerGMRF(ksp, z, &x, verbose_s);

        /* Share the sample */
        VecScatterBegin(scatter, x, x_scatter, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter, x, x_scatter, INSERT_VALUES, SCATTER_FORWARD);
        VecGetSubVector(x_scatter, is_separ, &x_separ);
        
        /* Correct */
        MatMult(Q_separ, x_separ, y);
        KSPSolve(ksp_exten, y, w);
        VecPointwiseMult(w, w, w);
        VecAXPY(d, 1.0/n_sample, w);
        VecRestoreSubVector(x_scatter, is_separ, &x_separ);
    }

    /* Assemble the solution */
    Vec d_origi, d_global;
    const int* is_array;
    double* d_array;
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, n, &d_global);
    VecGetSubVector(d, is_origi, &d_origi);
    ISGetIndices(is_local, &is_array);
    VecGetArray(d_origi, &d_array);
    for(int i=0; i<n_local; i++) VecSetValue(d_global, is_array[i], d_array[i], INSERT_VALUES);
    VecAssemblyBegin(d_global);
    VecAssemblyEnd(d_global);
    ISRestoreIndices(is_origi, &is_array);
    VecRestoreArray(d_origi, &d_array);

    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/out", FILE_MODE_WRITE, &viewer);
    VecView(d_global, viewer);

                /* Memory leak somewhere here... scatter? subvector? */
            double log;
            PetscMallocGetCurrentUsage(&log);
            if(!rank) printf("----------------- %f \n", log);


    /* Clean up */
    MatDestroy(&K1);
    MatDestroy(&K2);
    MatDestroy(&K3);
    MatDestroy(&J0);
    MatDestroy(&J1);
    MatDestroy(&J2);
    // MatDestroy(&J0_local);
    // MatDestroy(&J1_local);
    // MatDestroy(&J2_local);
    MatDestroy(&J0_global);
    MatDestroy(&J1_global);
    MatDestroy(&J2_global);
    MatDestroy(&J0_exten);
    MatDestroy(&J1_exten);
    MatDestroy(&J2_exten);
    MatDestroy(&J0_separ);
    MatDestroy(&J1_separ);
    MatDestroy(&J2_separ);
    // MatDestroy(&Q_local);
    // MatDestroy(&Q1_local);
    // MatDestroy(&Q2_local);
    MatDestroy(&Q);
    MatDestroy(&Q1);
    MatDestroy(&Q2);
    MatDestroy(&Q_exten);
    MatDestroy(&Q1_exten);
    MatDestroy(&Q2_exten);
    MatDestroy(&Q_separ);
    MatDestroy(&Q1_separ);
    MatDestroy(&Q2_separ);
    MatDestroy(&L);
    MatDestroy(&II);
    MatDestroy(&S);

    VecDestroy(&d);
    VecDestroy(&d_origi);
    VecDestroy(&d_global);
    VecDestroy(&x);
    VecDestroy(&x_scatter);
    // VecDestroy(&x_separ);
    VecDestroy(&z);
    VecDestroy(&y);
    VecDestroy(&w);

    ISDestroy(&ist_local);
    ISDestroy(&ist_ghost);
    ISDestroy(&ist_exten);
    ISDestroy(&ist_separ);
    ISDestroy(&is_local);
    ISDestroy(&is_ghost);
    ISDestroy(&is_exten);
    ISDestroy(&is_separ);
    ISDestroy(&is_direct);
    ISDestroy(&is_origi);

    KSPDestroy(&ksp);
    KSPDestroy(&ksp_exten);
    PetscViewerDestroy(&viewer);
    VecScatterDestroy(&scatter);
    InvGraphDestroy(graph, verbose);


    /* Finalize */
    PetscFinalize();
    MPI_Finalize();


    return 0;
}