#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <petsc.h>
#include <petscmat.h>
#include <mpi.h>

#include "inv_matrix.h"
#include "inv_graph.h"
#include "inv_shell.h"
#include "inv_sampler.h"
#include "inv_random.h"


int main(int argc, char **argv){


    /* Set parameters from the options */
    int max_niter = 1000, n_sample = 100, n_neighbor = 1, selected = 1, verbose = 0, verbose_s = 0, profile = 0;
    double tau_y = 1e-5, tau_b = 1e-5;
    double t_start = 0, t_end = 0, mem_start = 0, mem_end = 0;
    for(int i=0; i<argc; i++){
        if(!strcmp(argv[i], "-ns")){
            n_sample = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nmax")){
            max_niter = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nn")){
            n_neighbor = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-sel")){
            selected = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-tauy")){
            tau_y = atof(argv[i+1]);
        } else if(!strcmp(argv[i], "-taub")){
            tau_b = atof(argv[i+1]);
        } else if(!strcmp(argv[i], "-v")){
            verbose = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-vs")){
            verbose_s = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-prof")){
            profile = atoi(argv[i+1]);
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
    profile = (!rank) && profile;

    /* ---------------------------------------------------------------- */
    /* -------------------- PREPARATORY PHASE ------------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nPREPARATORY PHASE\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Load the matrices */
    if(verbose) printf("\tLoading the matrices...\n");
    PetscViewer viewer;
    Mat JJ0, JJ1, JJ2, J0, J1, J2, K1, K2, K3;
    MatCreate(PETSC_COMM_WORLD, &JJ0);
    MatCreate(PETSC_COMM_WORLD, &JJ1);
    MatCreate(PETSC_COMM_WORLD, &JJ2);
    MatCreate(PETSC_COMM_SELF, &J0);
    MatCreate(PETSC_COMM_SELF, &J1);
    MatCreate(PETSC_COMM_SELF, &J2);
    MatCreate(PETSC_COMM_SELF, &K1);
    MatCreate(PETSC_COMM_SELF, &K2);
    MatCreate(PETSC_COMM_SELF, &K3);
    MatSetType(JJ0, MATMPIAIJ);
    MatSetType(JJ1, MATMPIAIJ);
    MatSetType(JJ2, MATMPIAIJ);
    MatSetType(J0, MATSEQAIJ);
    MatSetType(J1, MATSEQAIJ);
    MatSetType(J2, MATSEQAIJ);
    MatSetType(K1, MATSEQAIJ);
    MatSetType(K2, MATSEQAIJ);
    MatSetType(K3, MATSEQAIJ);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/J0", FILE_MODE_READ, &viewer);
    MatLoad(JJ0, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/J1", FILE_MODE_READ, &viewer);
    MatLoad(JJ1, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/J2", FILE_MODE_READ, &viewer);
    MatLoad(JJ2, viewer);
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


    /* Load and partition the temporal graph */
    if(verbose) printf("\tLoading the graph...\n");
    int verbose_g = 0;
    verbose_g = verbose_g && verbose;
    InvGraph* graph = InvGraphCreate("data/J2", verbose_g);
    InvGraphPartition(graph, size, verbose_g);
    InvGraphNeighbor(graph, n_neighbor, verbose_g);
    // InvGraphPrint(graph, verbose_g);


    /* Set dimensions */
    if(verbose) printf("\tSetting dimensions...\n");
    int ns, nt, n;
    int nt_local, nt_ghost, nt_exten, nt_separ;
    int n_local, n_exten;
    MatGetSize(K3, &ns, &ns);
    nt_local = graph->offset[1] - graph->offset[0];
    nt_ghost = graph->offset[2] - graph->offset[1];
    nt_exten = graph->offset[2] - graph->offset[0];
    nt_separ = graph->offset[3] - graph->offset[2];
    nt = graph->n_vert;
    n_local = nt_local * ns;
    n_exten = nt_exten * ns;
    n = ns * nt;
    if(verbose) printf("\tProblem size: ns=%i, nt=%i, n=%i\n", ns, nt, n);
    if(verbose) printf("\tPartition sizes: nt_local=%i, nt_ghost=%i, nt_exten=%i, nt_separ=%i\n", nt_local, nt_ghost, nt_exten, nt_separ);
    if(verbose) printf("\tWorkload per process: n_local=%i, n_exten=%i\n", n_local, n_exten);


    /* Create index sets */
    if(verbose) printf("\tCreating index sets...\n");
    IS ist_local, ist_ghost, ist_exten, ist_separ, is_local, is_ghost, is_exten, is_separ, is_local2;
    ISCreateGeneral(PETSC_COMM_SELF, nt_local, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &ist_local);
    ISCreateGeneral(PETSC_COMM_SELF, nt_ghost, graph->vert + graph->offset[1], PETSC_COPY_VALUES, &ist_ghost);
    ISCreateGeneral(PETSC_COMM_SELF, nt_exten, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &ist_exten);
    ISCreateGeneral(PETSC_COMM_SELF, nt_separ, graph->vert + graph->offset[2], PETSC_COPY_VALUES, &ist_separ);
    ISCreateBlock(PETSC_COMM_SELF, ns, nt_local, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &is_local);
    ISCreateBlock(PETSC_COMM_SELF, ns, nt_ghost, graph->vert + graph->offset[1], PETSC_COPY_VALUES, &is_ghost);
    ISCreateBlock(PETSC_COMM_SELF, ns, nt_exten, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &is_exten);
    ISCreateBlock(PETSC_COMM_SELF, ns, nt_separ, graph->vert + graph->offset[2], PETSC_COPY_VALUES, &is_separ);
    ISCreateStride(PETSC_COMM_SELF, n_local, 0, 1, &is_local2);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(verbose) printf("\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(verbose) printf("\tMemory allocated:\t%f bytes\n", mem_end - mem_start);


    /* ---------------------------------------------------------------- */
    /* -------------------- MATRIX ASSEMBLY PHASE --------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nMATRIX ASSEMBLY PHASE\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Assemble the global MPI matrix */
    if(verbose) printf("\tAssembling the global matrix...\n");
    int istart, iend;
    Mat QQ, QQ1, QQ2;
    MatMPIAIJKron(JJ0, K3, &QQ);
    MatMPIAIJKron(JJ1, K2, &QQ1);
    MatMPIAIJKron(JJ2, K1, &QQ2);
    MatAXPY(QQ, 1.0, QQ1, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(QQ, 1.0, QQ2, DIFFERENT_NONZERO_PATTERN);
    MatGetOwnershipRange(QQ, &istart, &iend);
    for(int i=istart; i<iend; i++) MatSetValue(QQ, i, i, tau_y, ADD_VALUES);
    MatAssemblyBegin(QQ, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(QQ, MAT_FINAL_ASSEMBLY);


    /* Assemble the local extended matrix */
    if(verbose) printf("\tAssembling the local matrix...\n");
    Mat J0_exten, J1_exten, J2_exten, Q_exten, Q1_exten, Q2_exten;
    MatCreateSubMatrix(J0, ist_exten, ist_exten, MAT_INITIAL_MATRIX, &J0_exten);
    MatCreateSubMatrix(J1, ist_exten, ist_exten, MAT_INITIAL_MATRIX, &J1_exten);
    MatCreateSubMatrix(J2, ist_exten, ist_exten, MAT_INITIAL_MATRIX, &J2_exten);
    MatSeqAIJKron(J0_exten, K3, MAT_INITIAL_MATRIX, &Q_exten);
    MatSeqAIJKron(J1_exten, K2, MAT_INITIAL_MATRIX, &Q1_exten);
    MatSeqAIJKron(J2_exten, K1, MAT_INITIAL_MATRIX, &Q2_exten);
    MatAXPY(Q_exten, 1.0, Q1_exten, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Q_exten, 1.0, Q2_exten, DIFFERENT_NONZERO_PATTERN);
    MatGetOwnershipRange(Q_exten, &istart, &iend);
    for(int i=istart; i<iend; i++) MatSetValue(Q_exten, i, i, tau_y, ADD_VALUES);
    MatAssemblyBegin(Q_exten, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q_exten, MAT_FINAL_ASSEMBLY);


    /* Assemble the separator portion of the matrix */
    if(verbose) printf("\tAssembling the separator portion of the matrix...\n");
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


    /* Assemble the precision matrix with an intercept block */
    if(verbose) printf("\tAssembling the matrix with covariates...\n");
    Mat QQQ, QQ_ub, QQ_bu, QQ_bb;
    MatConcatenateIntercept(QQ, &QQQ, tau_y, tau_b);
    MatNestGetSubMat(QQQ, 0, 1, &QQ_ub);
    MatNestGetSubMat(QQQ, 1, 0, &QQ_bu);
    MatNestGetSubMat(QQQ, 1, 1, &QQ_bb);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(verbose) printf("\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(verbose) printf("\tMemory allocated:\t%f bytes\n", mem_end - mem_start);


    /* ---------------------------------------------------------------- */
    /* -------------------- DIRECT SOLVE PHASE ------------------------ */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nDIRECT SOLVE PHASE\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Factor the extended matrix */
    if(verbose) printf("\tFactoring the extended matrix...\n");
    Vec d;
    Mat L;
    IS is_factor;
    VecCreateSeq(PETSC_COMM_SELF, n_exten, &d);
    MatSetOption(Q_exten, MAT_SYMMETRIC, PETSC_TRUE);
    MatGetOrdering(Q_exten, MATORDERINGNATURAL, &is_factor, &is_factor);
    MatGetFactor(Q_exten, MATSOLVERMUMPS, MAT_FACTOR_CHOLESKY, &L);
    MatCholeskyFactorSymbolic(L, Q_exten, is_factor, NULL);
    MatCholeskyFactorNumeric(L, Q_exten, NULL);


    if(selected){
        /* Selected inversion */
        if(verbose) printf("\tSelected inversion...\n");
        Mat D;
        MatCreate(PETSC_COMM_SELF, &D);
        MatSetSizes(D, n_exten, n_exten, PETSC_DECIDE, PETSC_DECIDE);
        MatSetType(D, MATSEQAIJ);
        MatSetFromOptions(D);
        MatSetUp(D);
        for(int i=0; i<n_exten; i++) MatSetValue(D, i, i, 1.0, INSERT_VALUES);
        MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY);
        MatMumpsGetInverseTranspose(L, D);
        MatGetDiagonal(D, d);
        MatDestroy(&D);
    } else {
        /* Dense inversion */
        if(verbose) printf("\tDense inversion...\n");
        Mat II, S;
        MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, n_exten, n_exten, NULL, &II);
        MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, n_exten, n_exten, NULL, &S);
        for(int i=0; i<n_exten; i++) MatSetValue(II, i, i, 1.0, INSERT_VALUES);
        MatAssemblyBegin(II, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(II, MAT_FINAL_ASSEMBLY);
        MatMatSolve(L, II, S);
        MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
        MatGetDiagonal(S, d);
        MatDestroy(&II);
        MatDestroy(&S);
    }


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(verbose) printf("\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(verbose) printf("\tMemory allocated:\t%f bytes\n", mem_end - mem_start);


    /* ---------------------------------------------------------------- */
    /* -------------------- SAMPLING CORRECTION PHASE ----------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nSAMPLING CORRECTION PHASE\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Create a solver for sampling */
    if(verbose) printf("\tSetting up a sampling solver...\n");
    InvShellPC* shell;
    KSP ksp_sampling;
    PC pc_sampling;
    KSPCreate(PETSC_COMM_WORLD, &ksp_sampling);
    KSPSetOperators(ksp_sampling, QQ, QQ);
    KSPSetTolerances(ksp_sampling, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, max_niter);
    KSPSetComputeEigenvalues(ksp_sampling, PETSC_TRUE);
    KSPSetPCSide(ksp_sampling, PC_SYMMETRIC);
    KSPGetPC(ksp_sampling, &pc_sampling);

    if(verbose) printf("\tSetting up a shell preconditioner...\n");
    PCSetType(pc_sampling, PCSHELL);
    InvShellCreate(&shell);
    PCShellSetApply(pc_sampling, InvShellApply);
    PCShellSetApplyTranspose(pc_sampling, InvShellApplyTranspose);
    PCShellSetApplyBA(pc_sampling, InvShellApplyBA);
    PCShellSetApplySymmetricLeft(pc_sampling, InvShellApplyLeft);
    PCShellSetApplySymmetricRight(pc_sampling, InvShellApplyRight);
    PCShellSetContext(pc_sampling, shell);
    PCShellSetDestroy(pc_sampling, InvShellDestroy);
    PCShellSetName(pc_sampling, "shell");
    InvShellSetup(pc_sampling, QQ, PETSC_COMM_WORLD);


    /* Create a solver for correction */
    if(verbose) printf("\tSetting up a correction solver...\n");
    KSP ksp_correction;
    PC pc_correction;
    KSPCreate(PETSC_COMM_SELF, &ksp_correction);
    KSPSetOperators(ksp_correction, Q_exten, Q_exten);
    KSPSetType(ksp_correction, KSPGMRES);
    KSPGetPC(ksp_correction, &pc_correction);
    PCSetType(pc_correction, PCBJACOBI);
    PCSetFromOptions(pc_correction);
    KSPSetTolerances(ksp_correction, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, max_niter);
    KSPSetFromOptions(ksp_correction);


    /* Prepare sampling vectors */
    if(verbose) printf("\tSetting up vectors...\n");
    pcg64_random_t rng = InvRandomCreate(0, 0, 1);
    VecScatter scatter;
    Vec x, x_scatter, x_separ, z, v, w;
    int n_local_global;
    MatGetLocalSize(QQ, &n_local_global, &n_local_global);
    VecCreateMPI(PETSC_COMM_WORLD, n_local_global, PETSC_DETERMINE, &x);
    VecDuplicate(x, &z);
    VecScatterCreateToAll(x, &scatter, &x_scatter);
    VecCreateSeq(PETSC_COMM_SELF, n_exten, &v);
    VecDuplicate(v, &w);


    /* Sampling correction */
    if(verbose) printf("\tSampling correction...\n");
    for(int i=0; i<n_sample; i++){

        /* Sampling */
        InvSamplerStdNormal(&rng, &z, verbose_s);
        InvSamplerGMRF(ksp_sampling, z, &x, verbose_s);

        /* Share the sample */
        VecScatterBegin(scatter, x, x_scatter, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter, x, x_scatter, INSERT_VALUES, SCATTER_FORWARD);
        VecGetSubVector(x_scatter, is_separ, &x_separ);
        
        /* Correct */
        MatMult(Q_separ, x_separ, v);
        // MatSolve(L, y, w);
        KSPSolve(ksp_correction, v, w);  /* No need for solve with KSP, can use MatSolve(L) instead... */
        VecPointwiseMult(w, w, w);
        VecAXPY(d, 1.0/n_sample, w);
        VecRestoreSubVector(x_scatter, is_separ, &x_separ);
    }


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(verbose) printf("\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(verbose) printf("\tMemory allocated:\t%f bytes\n", mem_end - mem_start);


    /* ---------------------------------------------------------------- */
    /* -------------------- COVARIATES CORRECTION PHASE --------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nCOVARIATES CORRECTION PHASE\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Create a solver for covariates correction */
    if(verbose) printf("\tSetting up a solver for covariates correction...\n");
    int n_full, n_local_full;
    Vec e, s;
    KSP ksp_full;
    PC pc_full;
    MatGetLocalSize(QQQ, &n_local_full, &n_local_full);
    MatGetSize(QQQ, &n_full, &n_full);
    VecCreateMPI(PETSC_COMM_WORLD, n_local_full, PETSC_DETERMINE, &e);
    VecDuplicate(e, &s);
    VecSet(e, 0);
    VecSetValue(e, n_full-1, 1.0, INSERT_VALUES);
    VecAssemblyBegin(e);
    VecAssemblyEnd(e);

    KSPCreate(PETSC_COMM_WORLD, &ksp_full);
    KSPSetOperators(ksp_full, QQQ, QQQ);
    KSPGetPC(ksp_full, &pc_full);
    PCSetType(pc_full, PCFIELDSPLIT);
    PCFieldSplitSetType(pc_full, PC_COMPOSITE_ADDITIVE);
    PCFieldSplitSetSchurFactType(pc_full, PC_FIELDSPLIT_SCHUR_FACT_DIAG);
    KSPSetTolerances(ksp_full, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetUp(ksp_full);


    /* Solve and compute correction */
    if(verbose) printf("\tSolving and computing the correction...\n");
    double sigma_intercept;
    KSPSolve(ksp_full, e, s);
    VecPointwiseMult(e, s, e);
    VecPointwiseMult(s, s, s);
    VecSum(e, &sigma_intercept);
    VecScale(s, 1.0/sigma_intercept);


    /* Assemble the estimate conditioned on the intercept */
    if(verbose) printf("\tAssembling the approximation for the diagonal...\n");
    Vec d_local, d_global;
    const int* is_array;
    double* d_array;
    VecCreateMPI(PETSC_COMM_WORLD, n_local_full, PETSC_DETERMINE, &d_global);
    VecGetSubVector(d, is_local2, &d_local);
    ISGetIndices(is_local, &is_array);
    VecGetArray(d_local, &d_array);
    for(int i=0; i<n_local; i++) VecSetValue(d_global, is_array[i], d_array[i], INSERT_VALUES);
    VecAssemblyBegin(d_global);
    VecAssemblyEnd(d_global);
    ISRestoreIndices(is_local2, &is_array);
    VecRestoreArray(d_local, &d_array);
    VecAXPY(d_global, 1.0, s);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(verbose) printf("\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(verbose) printf("\tMemory allocated:\t%f bytes\n", mem_end - mem_start);


    /* ---------------------------------------------------------------- */
    /* -------------------- SOLVING FOR THE MEAN ---------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nSOLVING FOR THE MEAN\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Load the data */
    if(verbose) printf("\tLoading the data...\n");
    Vec y, b, mu;
    double* y_array;
    double y_sum;
    VecDuplicate(s, &b);
    VecDuplicate(s, &mu);
    VecCreate(PETSC_COMM_SELF, &y);
    if(!rank){
        PetscViewerBinaryOpen(PETSC_COMM_SELF, "data/y", FILE_MODE_READ, &viewer);
        VecLoad(y, viewer);
        VecSum(y, &y_sum);
        VecGetArray(y, &y_array);
        for(int i=0; i<n; i++) VecSetValue(b, i, tau_y*y_array[i], INSERT_VALUES);
        VecSetValue(b, n_full-1, tau_y*y_sum, INSERT_VALUES);
        VecRestoreArray(y, &y_array);
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);


    /* Solve for the mean */
    if(verbose) printf("\tSolving a system...\n");
    KSPSolve(ksp_full, b, mu);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(verbose) printf("\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(verbose) printf("\tMemory allocated:\t%f bytes\n", mem_end - mem_start);


    /* ---------------------------------------------------------------- */
    /* -------------------- FINALIZATION PHASE ------------------------ */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nFINALIZATION PHASE\n");


    /* Save the output */
    if(verbose) printf("\tSaving results...\n");
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/d", FILE_MODE_WRITE, &viewer);
    VecView(d_global, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/mu", FILE_MODE_WRITE, &viewer);
    VecView(mu, viewer);


    /* Clean up */
    if(verbose) printf("\tCleaning up...\n");
    MatDestroy(&K1);
    MatDestroy(&K2);
    MatDestroy(&K3);
    MatDestroy(&J0);
    MatDestroy(&J1);
    MatDestroy(&J2);
    MatDestroy(&JJ0);
    MatDestroy(&JJ1);
    MatDestroy(&JJ2);
    MatDestroy(&J0_exten);
    MatDestroy(&J1_exten);
    MatDestroy(&J2_exten);
    MatDestroy(&J0_separ);
    MatDestroy(&J1_separ);
    MatDestroy(&J2_separ);
    MatDestroy(&QQQ);
    MatDestroy(&QQ);
    MatDestroy(&QQ1);
    MatDestroy(&QQ2);
    MatDestroy(&Q_exten);
    MatDestroy(&Q1_exten);
    MatDestroy(&Q2_exten);
    MatDestroy(&Q_separ);
    MatDestroy(&Q1_separ);
    MatDestroy(&Q2_separ);
    MatDestroy(&QQ_ub);
    MatDestroy(&QQ_bu);
    MatDestroy(&QQ_bb);
    MatDestroy(&L);

    VecDestroy(&d);
    VecDestroy(&d_local);
    VecDestroy(&d_global);
    VecDestroy(&x);
    VecDestroy(&x_scatter);
    // VecDestroy(&x_separ);
    VecDestroy(&z);
    VecDestroy(&v);
    VecDestroy(&w);
    VecDestroy(&e);
    VecDestroy(&s);
    VecDestroy(&y);
    VecDestroy(&b);
    VecDestroy(&mu);

    ISDestroy(&ist_local);
    ISDestroy(&ist_ghost);
    ISDestroy(&ist_exten);
    ISDestroy(&ist_separ);
    ISDestroy(&is_local);
    ISDestroy(&is_ghost);
    ISDestroy(&is_exten);
    ISDestroy(&is_separ);
    ISDestroy(&is_local2);
    ISDestroy(&is_factor);

    KSPDestroy(&ksp_sampling);
    KSPDestroy(&ksp_correction);
    KSPDestroy(&ksp_full);

    PetscViewerDestroy(&viewer);
    VecScatterDestroy(&scatter);
    InvGraphDestroy(graph, verbose_g);


    /* Finalize */
    if(verbose) printf("\tFinalizing...\n");
    PetscFinalize();
    MPI_Finalize();


    return 0;
}