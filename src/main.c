#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <petsc.h>
#include <mpi.h>

#include "inv_is.h"
#include "inv_vector.h"
#include "inv_matrix.h"
#include "inv_random.h"
#include "inv_sampler.h"
#include "inv_shell.h"


int main(int argc, char **argv){


    /* Set parameters from the options */
    int n_iter = 1000, n_sample = 100, n_neighbor = 1, verbose = 0, profile = 0;
    double tau_y = 1e-2, tau_b = 1e-5;
    for(int i=0; i<argc; i++){
        if(!strcmp(argv[i], "-ns")){
            n_sample = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-ni")){
            n_iter = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nn")){
            n_neighbor = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-v")){
            verbose = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-p")){
            profile = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-tauy")){
            tau_y = atof(argv[i+1]);
        } else if(!strcmp(argv[i], "-taub")){
            tau_b = atof(argv[i+1]);
        }
        argv[i] = 0;
    }


    /* Initialize MPI and PETSc */
    int rank, size, last;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    PetscInitialize(&argc, &argv, NULL, NULL);
    verbose = (!rank) && verbose;
    profile = verbose && profile ;
    last = (rank==size-1);
    double t_init = 0, t_start = 0, t_end = 0, mem_init = 0, mem_start = 0, mem_end = 0;


    /* ---------------------------------------------------------------- */
    /* ---------------------------- INPUT ----------------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nINPUT\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Set global dimensions */
    if(verbose) printf("\tSetting global dimensions...\n");
    int nb;
    int nt, ns, nu;
    int mt, ms, mu;
    MatGetSizeLoad("data/At", &mt, &nt);
    MatGetSizeLoad("data/As", &ms, &ns);
    MatGetSizeLoad("data/Ab", &mu, &nb);
    nu = nt * ns;
    mu = mt * ms;


    /* Load structure matrices */
    if(verbose) printf("\tLoading matrices...\n");
    Mat J0, J1, J2, At;
    Mat K1, K2, K3, As;
    MatCreateLoad(PETSC_COMM_WORLD, MATMPIAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/J0", &J0);
    MatCreateLoad(PETSC_COMM_WORLD, MATMPIAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/J1", &J1);
    MatCreateLoad(PETSC_COMM_WORLD, MATMPIAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/J2", &J2);
    MatCreateLoad(PETSC_COMM_WORLD, MATMPIAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/At", &At);
    MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/K1", &K1);
    MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/K2", &K2);
    MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/K3", &K3);
    MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/As", &As);


    /* Set local dimensions */
    if(verbose) printf("\tSetting local dimensions...\n");
    int nt_local, nu_local;
    int mt_local, mu_local;
    MatGetLocalSize(At, &mt_local, &nt_local);
    nu_local = nt_local * ns;
    mu_local = mt_local * ms;


    /* Load the data */
    if(verbose) printf("\tLoading the data...\n");
    Mat Ab;
    Vec y;
    IS is_na;
    MatCreateLoad(PETSC_COMM_WORLD, MATMPIDENSE, mu_local, nb*last, PETSC_DETERMINE, nb, "data/Ab", &Ab);
    VecCreateLoad(PETSC_COMM_WORLD, VECMPI, mu_local, PETSC_DETERMINE, "data/y", &y);
    ISCreateLoad(PETSC_COMM_SELF, "data/isna", &is_na);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- PARTITIONING ---------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nPARTITIONING\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Compute the ghost and separator nodes */
    if(verbose) printf("\tComputing the ghost and separator nodes...\n");
    IS isnt_local, isnt_ghost, isnt_exten, isnt_separ;
    MatGetOwnershipISRows(PETSC_COMM_WORLD, J2, &isnt_local);
    MatGetOwnershipISRows(PETSC_COMM_WORLD, J2, &isnt_exten);
    MatIncreaseOverlap(J2, 1, &isnt_exten, n_neighbor);
    ISDuplicate(isnt_exten, &isnt_ghost);
    ISDuplicate(isnt_exten, &isnt_separ);
    ISCopy(isnt_exten, isnt_ghost);
    ISCopy(isnt_exten, isnt_separ);
    MatIncreaseOverlap(J2, 1, &isnt_separ, 1);
    ISDifference(isnt_exten, isnt_local, &isnt_ghost);
    ISDifference(isnt_separ, isnt_exten, &isnt_separ);

    int nt_ghost, nt_exten, nt_separ;
    int mt_ghost, mt_exten, mt_separ;
    int nu_ghost, nu_exten, nu_separ;
    int mu_ghost, mu_exten, mu_separ;
    ISGetLocalSize(isnt_ghost, &nt_ghost);
    ISGetLocalSize(isnt_exten, &nt_exten);
    ISGetLocalSize(isnt_separ, &nt_separ);
    mt_ghost = nt_ghost;    //  might change if At is not diagonal
    mt_exten = nt_exten;    //  ...
    mt_separ = nt_separ;    //  ...
    nu_ghost = nt_ghost * ns;
    nu_exten = nt_exten * ns;
    nu_separ = nt_separ * ns;
    mu_ghost = mt_ghost * ms;
    mu_exten = mt_exten * ms;
    mu_separ = mt_separ * ms;

    
    /* Compute the ghost and separator blocks */
    if(verbose) printf("\tComputing the ghost and separator blocks...\n");
    IS isnu_local, isnu_ghost, isnu_exten, isnu_separ;
    IS ismu_local, ismu_ghost, ismu_exten, ismu_separ;
    ISCreateBlockIS(PETSC_COMM_WORLD, isnt_local, ns, &isnu_local);
    ISCreateBlockIS(PETSC_COMM_WORLD, isnt_ghost, ns, &isnu_ghost);
    ISCreateBlockIS(PETSC_COMM_WORLD, isnt_exten, ns, &isnu_exten);
    ISCreateBlockIS(PETSC_COMM_WORLD, isnt_separ, ns, &isnu_separ);
    ISCreateBlockIS(PETSC_COMM_WORLD, isnt_local, ms, &ismu_local);
    ISCreateBlockIS(PETSC_COMM_WORLD, isnt_ghost, ms, &ismu_ghost);
    ISCreateBlockIS(PETSC_COMM_WORLD, isnt_exten, ms, &ismu_exten);
    ISCreateBlockIS(PETSC_COMM_WORLD, isnt_separ, ms, &ismu_separ);
    ISToGeneral(isnu_local);
    ISToGeneral(isnu_ghost);
    ISToGeneral(isnu_exten);
    ISToGeneral(isnu_separ);
    ISToGeneral(ismu_local);
    ISToGeneral(ismu_ghost);
    ISToGeneral(ismu_exten);
    ISToGeneral(ismu_separ);


    /* Map NA indices to ghost and separator blocks */
    if(verbose) printf("\tMapping NA observations...\n");
    IS ismu_local_na, ismu_ghost_na, ismu_exten_na, ismu_separ_na;
    ISLocalToGlobalMapping ltog_ismu_local, ltog_ismu_ghost, ltog_ismu_exten, ltog_ismu_separ;    
    ISLocalToGlobalMappingCreateIS(ismu_local, &ltog_ismu_local);
    ISLocalToGlobalMappingCreateIS(ismu_ghost, &ltog_ismu_ghost);
    ISLocalToGlobalMappingCreateIS(ismu_exten, &ltog_ismu_exten);
    ISLocalToGlobalMappingCreateIS(ismu_separ, &ltog_ismu_separ);
    ISLocalToGlobalMappingSetType(ltog_ismu_local, ISLOCALTOGLOBALMAPPINGHASH);
    ISLocalToGlobalMappingSetType(ltog_ismu_ghost, ISLOCALTOGLOBALMAPPINGHASH);
    ISLocalToGlobalMappingSetType(ltog_ismu_exten, ISLOCALTOGLOBALMAPPINGHASH);
    ISLocalToGlobalMappingSetType(ltog_ismu_separ, ISLOCALTOGLOBALMAPPINGHASH);
    ISGlobalToLocalMappingApplyIS(ltog_ismu_local, IS_GTOLM_DROP, is_na, &ismu_local_na);
    ISGlobalToLocalMappingApplyIS(ltog_ismu_ghost, IS_GTOLM_DROP, is_na, &ismu_ghost_na);
    ISGlobalToLocalMappingApplyIS(ltog_ismu_exten, IS_GTOLM_DROP, is_na, &ismu_exten_na);
    ISGlobalToLocalMappingApplyIS(ltog_ismu_separ, IS_GTOLM_DROP, is_na, &ismu_separ_na);


    /* Print problem dimensions */
    if(verbose) printf("\n\tGlobal dimensions:\tnt=%i, ns=%i, nu=%i, mt=%i, ms=%i, mu=%i\n", nt, ns, nu, mt, ms, mu);
    if(verbose) printf("\tLocal dimensions:\tnt=%i, ns=%i, nu=%i, mt=%i, ms=%i, mu=%i\n", nt_local, ns, nu_local, mt_local, ms, mu_local);
    if(verbose) printf("\tGhost dimensions:\tnt=%i, ns=%i, nu=%i, mt=%i, ms=%i, mu=%i\n", nt_ghost, ns, nu_ghost, mt_ghost, ms, mu_ghost);
    if(verbose) printf("\tExtended dimensions:\tnt=%i, ns=%i, nu=%i, mt=%i, ms=%i, mu=%i\n", nt_exten, ns, nu_exten, mt_exten, ms, mu_exten);
    if(verbose) printf("\tSeparator dimensions:\tnt=%i, ns=%i, nu=%i, mt=%i, ms=%i, mu=%i\n", nt_separ, ns, nu_separ, mt_separ, ms, mu_separ);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- ASSEMBLY -------------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nASSEMBLY\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Compute the global precision matrix */
    if(verbose) printf("\tComputing the global precision matrix...\n");
    Mat Au;
    Mat Quu, Quu1, Quu2, Quu3;
    Mat Qub;
    Mat Qbb, Qbb1;
    MatMPIAIJKron(J0, K3, &Quu);
    MatMPIAIJKron(J1, K2, &Quu1);
    MatMPIAIJKron(J2, K1, &Quu2);
    MatMPIAIJKron(At, As, &Au);
    MatZeroRowsIS(Au, is_na, 0.0, NULL, NULL);
    MatZeroRowsIS(Ab, is_na, 0.0, NULL, NULL);
    MatCreateDense(PETSC_COMM_WORLD, nb*last, nb*last, nb, nb, NULL, &Qbb);
    MatShift(Qbb, tau_b);
    MatTransposeMatMult(Au, Au, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Quu3);
    MatTransposeMatMult(Ab, Ab, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Qbb1);
    MatTransposeMatMult(Au, Ab, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Qub);
    MatAXPY(Quu, 1.0, Quu1, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Quu, 1.0, Quu2, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Quu, 1.0, Quu3, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Qbb, 1.0, Qbb1, SAME_NONZERO_PATTERN);


    /* Compute the local extended precision matrix */
    if(verbose) printf("\tComputing the local extended precision matrix...\n");
    Mat J0_exten, J1_exten, J2_exten, At_exten;
    MatCreateSubMatrixSeq(J0, isnt_exten, isnt_exten, &J0_exten);
    MatCreateSubMatrixSeq(J1, isnt_exten, isnt_exten, &J1_exten);
    MatCreateSubMatrixSeq(J2, isnt_exten, isnt_exten, &J2_exten);
    MatCreateSubMatrixSeq(At, isnt_exten, isnt_exten, &At_exten);   // if At is not diagonal, then extract At[,ismt_exten] (costly?)
                                                                    //  Ideally, keep At distributed according to ismt_local, then At^T * At will automatically collect the needed blocks
    Mat Au_exten;
    Mat Quu_exten, Quu1_exten, Quu2_exten, Quu3_exten;
    MatSeqAIJKron(J0_exten, K3, MAT_INITIAL_MATRIX, &Quu_exten);
    MatSeqAIJKron(J1_exten, K2, MAT_INITIAL_MATRIX, &Quu1_exten);
    MatSeqAIJKron(J2_exten, K1, MAT_INITIAL_MATRIX, &Quu2_exten);
    MatSeqAIJKron(At_exten, As, MAT_INITIAL_MATRIX, &Au_exten);
    MatZeroRowsIS(Au_exten, ismu_exten_na, 0.0, NULL, NULL);        //  if At is not diagonal, can set At[is_na,] to zero
    MatTransposeMatMult(Au_exten, Au_exten, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Quu3_exten);
    MatAXPY(Quu_exten, 1.0, Quu1_exten, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Quu_exten, 1.0, Quu2_exten, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Quu_exten, 1.0, Quu3_exten, DIFFERENT_NONZERO_PATTERN);


    /* Compute the separator extended precision matrix */
    if(verbose) printf("\tComputing the separator extended precision matrix...\n");
    Mat J0_separ, J1_separ, J2_separ, At_separ;
    MatCreateSubMatrixSeq(J0, isnt_exten, isnt_separ, &J0_separ);
    MatCreateSubMatrixSeq(J1, isnt_exten, isnt_separ, &J1_separ);
    MatCreateSubMatrixSeq(J2, isnt_exten, isnt_separ, &J2_separ);
    MatCreateSubMatrixSeq(At, isnt_exten, isnt_separ, &At_separ);   // if At is not diagonal, then extract At[,ismt_separ]...

    Mat Au_separ;
    Mat Quu_separ, Quu1_separ, Quu2_separ, Quu3_separ;
    MatSeqAIJKron(J0_separ, K3, MAT_INITIAL_MATRIX, &Quu_separ);
    MatSeqAIJKron(J1_separ, K2, MAT_INITIAL_MATRIX, &Quu1_separ);
    MatSeqAIJKron(J2_separ, K1, MAT_INITIAL_MATRIX, &Quu2_separ);
    MatSeqAIJKron(At_separ, As, MAT_INITIAL_MATRIX, &Au_separ);
    MatZeroRowsIS(Au_separ, ismu_separ_na, 0.0, NULL, NULL);        //  if At is not diagonal, can set At[is_na,] to zero
    MatTransposeMatMult(Au_exten, Au_separ, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Quu3_separ);
    MatAXPY(Quu_separ, 1.0, Quu1_separ, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Quu_separ, 1.0, Quu2_separ, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(Quu_separ, 1.0, Quu3_separ, DIFFERENT_NONZERO_PATTERN);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));

    
    /* ---------------------------------------------------------------- */
    /* ---------------------------- DIRECT SOLVE ---------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nDIRECT SOLVE\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Factor the extended matrix */
    if(verbose) printf("\tFactoring the extended matrix...\n");
    Vec d_exten;
    Mat L_exten;
    IS is_factor;
    VecCreateSeq(PETSC_COMM_SELF, nu_exten, &d_exten);
    MatSetOption(Quu_exten, MAT_SYMMETRIC, PETSC_TRUE);
    MatGetOrdering(Quu_exten, MATORDERINGNATURAL, &is_factor, &is_factor);
    MatGetFactor(Quu_exten, MATSOLVERMUMPS, MAT_FACTOR_CHOLESKY, &L_exten);
    MatCholeskyFactorSymbolic(L_exten, Quu_exten, is_factor, NULL);
    MatCholeskyFactorNumeric(L_exten, Quu_exten, NULL);


    /* Selected inversion */
    if(verbose) printf("\tSelected inversion...\n");
    Mat D_exten;
    MatCreate(PETSC_COMM_SELF, &D_exten);
    MatSetSizes(D_exten, nu_exten, nu_exten, PETSC_DECIDE, PETSC_DECIDE);
    MatSetType(D_exten, MATSEQAIJ);
    MatSetFromOptions(D_exten);
    MatSetUp(D_exten);
    for(int i=0; i<nu_exten; i++) MatSetValue(D_exten, i, i, 1.0, INSERT_VALUES);
    MatAssemblyBegin(D_exten, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(D_exten, MAT_FINAL_ASSEMBLY);
    MatMumpsGetInverseTranspose(L_exten, D_exten);
    MatGetDiagonal(D_exten, d_exten);
    MatDestroy(&D_exten);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- KSP SOLVERS ----------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nKSP SOLVERS\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Create a solver for sampling */
    if(verbose) printf("\tSetting up a sampling solver...\n");
    InvShellPC* shell;
    KSP ksp_sampling;
    PC pc_sampling;
    KSPCreate(PETSC_COMM_WORLD, &ksp_sampling);
    KSPSetOperators(ksp_sampling, Quu, Quu);
    KSPSetTolerances(ksp_sampling, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, n_iter);
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
    InvShellSetup(pc_sampling, Quu, PETSC_COMM_WORLD);


    /* Create a solver for correction */
    if(verbose) printf("\tSetting up a correction solver...\n");
    KSP ksp_correction;
    PC pc_correction;
    KSPCreate(PETSC_COMM_SELF, &ksp_correction);
    KSPSetOperators(ksp_correction, Quu_exten, Quu_exten);
    KSPSetType(ksp_correction, KSPGMRES);
    KSPGetPC(ksp_correction, &pc_correction);
    PCSetType(pc_correction, PCBJACOBI);
    PCSetFromOptions(pc_correction);
    KSPSetTolerances(ksp_correction, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, n_iter);
    KSPSetFromOptions(ksp_correction);


    /* Create a solver for covariates */
    KSP ksp_covariates;
    PC pc_covariates;
    KSPCreate(PETSC_COMM_WORLD, &ksp_covariates);
    KSPSetOperators(ksp_covariates, Quu, Quu);
    KSPSetType(ksp_covariates, KSPGMRES);
    KSPGetPC(ksp_covariates, &pc_covariates);
    PCSetType(pc_covariates, PCBJACOBI);
    PCSetFromOptions(pc_covariates);
    KSPSetTolerances(ksp_covariates, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, n_iter);
    KSPSetFromOptions(ksp_covariates);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- SAMPLING CORRECTION --------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nSAMPLING CORRECTION\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Prepare sampling vectors */
    if(verbose) printf("\tSetting up vectors...\n");
    pcg64_random_t rng = InvRandomCreate(0, 0, 1);
    VecScatter scatter;
    Vec x, x_separ, x_separ_seq, z, v_exten, w_exten;
    VecCreateMPI(PETSC_COMM_WORLD, nu_local, PETSC_DETERMINE, &x);
    VecCreateMPI(PETSC_COMM_WORLD, nu_separ, PETSC_DETERMINE, &x_separ);
    VecCreateSeq(PETSC_COMM_SELF, nu_exten, &v_exten);
    VecDuplicate(x, &z);
    VecDuplicate(v_exten, &w_exten);
    VecScatterCreate(x, isnu_separ, x_separ, NULL, &scatter);
    VecCreateLocalVector(x_separ, &x_separ_seq);
    

    /* Sampling correction */
    if(verbose) printf("\tSampling correction...\n");
    for(int i=0; i<n_sample; i++){

        /* Sampling */
        InvSamplerStdNormal(&rng, &z);
        InvSamplerGMRF(ksp_sampling, z, &x);

        /* Share the sample */
        VecScatterBegin(scatter, x, x_separ, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter, x, x_separ, INSERT_VALUES, SCATTER_FORWARD);
        VecGetLocalVector(x_separ, x_separ_seq);
        
        /* Correct */
        MatMult(Quu_separ, x_separ_seq, v_exten);
        // MatSolve(L_exten, v_exten, w_exten);     /* Can use MatSolve() instead of KSPSolve(), but MUMPS throws an error... */
        KSPSolve(ksp_correction, v_exten, w_exten);
        VecPointwiseMult(w_exten, w_exten, w_exten);
        VecAXPY(d_exten, 1.0/n_sample, w_exten);
        VecRestoreLocalVector(x_separ, x_separ_seq);
    }


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- COVARIATE CORRECTION -------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nCOVARIATE CORRECTION\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Compute the Schur's complement */
    if(verbose) printf("\tComputing the Schur's complement...\n");
    Mat S, C, D;
    MatDuplicate(Qub, MAT_DO_NOT_COPY_VALUES, &C);
    KSPMatSolve(ksp_covariates, Qub, C);
    MatTransposeMatMult(Qub, C, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &S);
    MatAYPX(S, -1.0, Qbb, SAME_NONZERO_PATTERN);
    MatDenseInvertLapack(S);
    MatMatMult(C, S, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &D);
    MatDensePointwiseMult(D, C);


    /* Compute the covariate correction */
    if(verbose) printf("\tComputing the covariates correction...\n");
    Vec r, s;
    VecDuplicate(x, &r);
    VecCreateMPI(PETSC_COMM_WORLD, nb*last, PETSC_DETERMINE, &s);
    MatGetRowSum(D, r);
    MatGetDiagonal(S, s);


    /* Assemble the diagonal approximation */
    if(verbose) printf("\tAssembling the diagonal approximation...\n");
    Vec d, d_local, d_concat[2], diag;
    IS isnu_local_local;
    ISLocalToGlobalMapping ltog_isnu_local;
    ISLocalToGlobalMappingCreateIS(isnu_local, &ltog_isnu_local);
    ISLocalToGlobalMappingSetType(ltog_isnu_local, ISLOCALTOGLOBALMAPPINGHASH);
    ISGlobalToLocalMappingApplyIS(ltog_isnu_local, IS_GTOLM_DROP, isnu_local, &isnu_local_local);
    VecGetSubVector(d_exten, isnu_local_local, &d_local);
    VecCreateSuperVector(d_local, &d);
    VecRestoreSubVector(d_exten, isnu_local_local, &d_local);
    VecAXPY(d, 1.0, r);
    d_concat[0] = d;
    d_concat[1] = s;
    VecConcatenate(2, d_concat, &diag, NULL);    


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- SOLVE FOR THE MEAN ---------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nSOLVE FOR THE MEAN\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Solve the system */
    if(verbose) printf("\tSolving the system...\n");
    Vec bu, bb, cu, cb, meanu, meanb, mean_concat[2], mean;
    MatZeroRowsIS(Ab, is_na, 0.0, x, y);
    VecDuplicate(d, &bu);
    VecDuplicate(s, &bb);
    VecDuplicate(d, &cu);
    VecDuplicate(s, &cb);
    VecDuplicate(d, &meanu);
    VecDuplicate(s, &meanb);
    MatMultTranspose(Au, y, bu);
    MatMultTranspose(Ab, y, bb);
    MatMultTranspose(Qub, bu, cb);
    VecAXPY(bb, -1.0, cb);
    MatMult(S, bb, meanb);
    MatMult(Qub, meanb, cu);
    VecAXPY(bu, -1.0, cu);
    KSPSolve(ksp_covariates, bu, meanu);
    mean_concat[0] = meanu;
    mean_concat[1] = meanb;
    VecConcatenate(2, mean_concat, &mean, NULL);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- CLEAN UP -------------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nCLEAN UP\n");
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Save the output */
    if(verbose) printf("\tSaving results...\n");
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/d", FILE_MODE_WRITE, &viewer);
    VecView(diag, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/mu", FILE_MODE_WRITE, &viewer);
    VecView(mean, viewer);


    /* Destroy matrices */
    if(verbose) printf("\tDestroying matrices...\n");
    MatDestroy(&J0);
    MatDestroy(&J1);
    MatDestroy(&J2);
    MatDestroy(&At);
    MatDestroy(&K1);
    MatDestroy(&K2);
    MatDestroy(&K3);
    MatDestroy(&As);
    MatDestroy(&Ab);

    MatDestroy(&Au);
    MatDestroy(&Quu);
    MatDestroy(&Quu1);
    MatDestroy(&Quu2);
    MatDestroy(&Quu3);
    MatDestroy(&Qub);
    MatDestroy(&Qbb);
    MatDestroy(&Qbb1);

    MatDestroy(&J0_exten);
    MatDestroy(&J1_exten);
    MatDestroy(&J2_exten);
    MatDestroy(&At_exten);

    MatDestroy(&Au_exten);
    MatDestroy(&Quu_exten);
    MatDestroy(&Quu1_exten);
    MatDestroy(&Quu2_exten);
    MatDestroy(&Quu3_exten);

    MatDestroy(&J0_separ);
    MatDestroy(&J1_separ);
    MatDestroy(&J2_separ);
    MatDestroy(&At_separ);

    MatDestroy(&Au_separ);
    MatDestroy(&Quu_separ);
    MatDestroy(&Quu1_separ);
    MatDestroy(&Quu2_separ);
    MatDestroy(&Quu3_separ);

    MatDestroy(&L_exten);
    MatDestroy(&D_exten);

    MatDestroy(&S);
    MatDestroy(&C);
    MatDestroy(&D);
    

    /* Destroy vectors */
    if(verbose) printf("\tDestroying vectors...\n");
    VecDestroy(&y);
    VecDestroy(&d_exten);

    VecDestroy(&x);
    VecDestroy(&x_separ);
    VecDestroy(&x_separ_seq);
    VecDestroy(&z);
    VecDestroy(&v_exten);
    VecDestroy(&w_exten);

    VecDestroy(&r);
    VecDestroy(&s);
    VecDestroy(&d);
    VecDestroy(&d_local);
    VecDestroy(&diag);

    VecDestroy(&bu);
    VecDestroy(&bb);
    VecDestroy(&cu);
    VecDestroy(&cb);
    VecDestroy(&meanu);
    VecDestroy(&meanb);
    VecDestroy(&mean);


    /* Destroy index sets */
    if(verbose) printf("\tDestroying index sets...\n");
    ISDestroy(&is_na);

    ISDestroy(&isnt_local);
    ISDestroy(&isnt_ghost);
    ISDestroy(&isnt_exten);
    ISDestroy(&isnt_separ);

    ISDestroy(&isnu_local);
    ISDestroy(&isnu_ghost);
    ISDestroy(&isnu_exten);
    ISDestroy(&isnu_separ);

    ISDestroy(&ismu_local);
    ISDestroy(&ismu_ghost);
    ISDestroy(&ismu_exten);
    ISDestroy(&ismu_separ);

    ISDestroy(&ismu_local_na);
    ISDestroy(&ismu_ghost_na);
    ISDestroy(&ismu_exten_na);
    ISDestroy(&ismu_separ_na);

    ISDestroy(&is_factor);
    ISDestroy(&isnu_local_local);

    ISLocalToGlobalMappingDestroy(&ltog_ismu_local);
    ISLocalToGlobalMappingDestroy(&ltog_ismu_ghost);
    ISLocalToGlobalMappingDestroy(&ltog_ismu_exten);
    ISLocalToGlobalMappingDestroy(&ltog_ismu_separ);
    ISLocalToGlobalMappingDestroy(&ltog_isnu_local);


    /* Destroy solvers and other objects */
    if(verbose) printf("\tDestroying solvers and other objects...\n");
    KSPDestroy(&ksp_sampling);
    KSPDestroy(&ksp_correction);
    KSPDestroy(&ksp_covariates);

    // InvShellDestroy(pc_sampling);
    // PCDestroy(&pc_correction);
    // PCDestroy(&pc_covariates);

    VecScatterDestroy(&scatter);
    PetscViewerDestroy(&viewer);


    /* Profiling checkpoint */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    // /* ---------------------------------------------------------------- */
    // /* ---------------------------- FINALIZE -------------------------- */
    // /* ---------------------------------------------------------------- */
    if(verbose) printf("\nFINALIZE\n");


    /* Final profiling */
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTotal time:\t\t%f sec\n", t_end - t_init);
    if(profile) printf("\tUnfreed memory:\t\t%i bytes\n", (int)(mem_end - mem_init));


    /* Finalize */
    if(verbose) printf("\tFinalizing...\n");
    PetscFinalize();
    MPI_Finalize();


    return 0;
}


    // /* ---------------------------------------------------------------- */
    // /* ---------------------------- INPUT ----------------------------- */
    // /* ---------------------------------------------------------------- */
    // if(verbose) printf("\nINPUT\n");
    // if(profile) PetscTime(&t_start);
    // if(profile) PetscMallocGetCurrentUsage(&mem_start);


    // /* Set problems dimensions */
    // if(verbose) printf("\tSetting dimensions...\n");
    // int nb;
    // int nt, ns, nu;
    // int mt, ms, mu;
    // MatGetSizeLoad("data/At", &mt, &nt);
    // MatGetSizeLoad("data/As", &ms, &ns);
    // MatGetSizeLoad("data/Ab", &mu, &nb);
    // nu = nt * ns;
    // mu = mt * ms;


    // /* Load and partition the temporal graph */
    // if(verbose) printf("\tLoading the data...\n");
    // Mat J0_pre, J1_pre, J2_pre, At_pre;
    // Mat J0, J1, J2, At;
    // Mat K1, K2, K3, As;
    // MatCreateLoad(PETSC_COMM_WORLD, MATMPIAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/J2", &J2);
    
    // if(verbose) printf("\tPartitioning the temporal graph...\n");
    // int part_count[size];
    // MatPartitioning part;
    // IS is_part, is_part_num, is_perm, is_perm_inv;
    // MatPartitioningCreate(PETSC_COMM_WORLD, &part);
    // MatPartitioningSetAdjacency(part, J2_pre);
    // MatPartitioningSetNParts(part, size);
    // MatPartitioningSetType(part, MATPARTITIONINGPARMETIS);
    // MatPartitioningApply(part, &is_part);
    // ISPartitioningCount(is_part, size, part_count);
    // ISPartitioningToNumbering(is_part, &is_part_num);
    // ISInvertPermutation(is_part_num, part_count[rank], &is_perm);
    // ISInvertPermutation(is_perm, PETSC_DECIDE, &is_perm_inv);
    
    // // // Mat J2;
    // // MatView(J2_pre, PETSC_VIEWER_STDOUT_WORLD);
    // // MatCreateSubMatrix(J2_pre, is_perm, is_perm, MAT_INITIAL_MATRIX, &J2);
    // // MatView(J2_pre, PETSC_VIEWER_STDOUT_WORLD);

    // // IS is_perm1;
    // // ISDuplicate(is_perm, &is_perm1);
    // // ISCopy(is_perm, is_perm1);
    // // MatIncreaseOverlap(J2, 1, &is_perm1, 1);
    // // ISView(is_perm, PETSC_VIEWER_STDOUT_WORLD);
    // // ISView(is_perm1, PETSC_VIEWER_STDOUT_WORLD);
    // // // if(!rank) for(int i=0; i<size; i++) printf("\t%i\n", part_count[i]);

    // IS isrow, iscol, isrow1;
    // int istart, iend;
    // int nlocal, mlocal;
    // MatGetOwnershipRange(J2_pre, &istart, &iend);
    // MatGetLocalSize(J2_pre, &nlocal, &mlocal);
    // ISCreateStride(PETSC_COMM_WORLD, nlocal, istart, 1, &isrow);
    // ISCreateStride(PETSC_COMM_WORLD, nlocal, istart, 1, &isrow1);
    // MatScale(J2_pre, 1/0.01989437);
    // for(int i=istart; i<iend; i++) MatSetValue(J2_pre, i, i, i/10.0, ADD_VALUES);
    // MatAssemblyBegin(J2_pre, MAT_FINAL_ASSEMBLY);
    // MatAssemblyEnd(J2_pre, MAT_FINAL_ASSEMBLY);
    // // ISCreateStride(PETSC_COMM_WORLD, ng)
    // ISView(isrow, PETSC_VIEWER_STDOUT_WORLD);
    // MatIncreaseOverlap(J2_pre, 1, &isrow, 1);
    // ISView(isrow, PETSC_VIEWER_STDOUT_WORLD);

    // MatCreateSubMatrix(J2_pre, isrow, isrow, MAT_INITIAL_MATRIX, &J2);
    // MatPermute(J2_pre, isrow, isrow, &J0);
    // IS isr[1], isc[0];  isr[0] = isrow; isc[0] = isrow;
    // Mat *JJ[1];
    // MatCreateSubMatrices(J2_pre, 1, &isrow, &isrow, MAT_INITIAL_MATRIX, JJ);
    // J1 = *JJ[0];
    // Mat JJJ;
    // MatCreateSubMatrixSeq(J2_pre, isrow, isrow, &JJJ);
    // MatView(J2_pre, PETSC_VIEWER_STDOUT_WORLD);
    // MatView(J2, PETSC_VIEWER_STDOUT_WORLD);
    // MatView(J0, PETSC_VIEWER_STDOUT_WORLD);
    // if(rank==2)MatView(J1, PETSC_VIEWER_STDOUT_SELF);
    // if(rank==2)MatView(JJJ, PETSC_VIEWER_STDOUT_SELF);

    // ISLocalToGlobalMapping itog;
    // ISLocalToGlobalMappingCreateIS(isrow1, &itog);
    // ISLocalToGlobalMappingSetType(itog, ISLOCALTOGLOBALMAPPINGHASH);
    // ISLocalToGlobalMappingView(itog, PETSC_VIEWER_STDOUT_WORLD);

    // ISLocalToGlobalMapping ltog;
    // ISLocalToGlobalMappingCreateIS(isrow, &ltog);
    // ISLocalToGlobalMappingSetType(ltog, ISLOCALTOGLOBALMAPPINGHASH);
    // ISLocalToGlobalMappingView(ltog, PETSC_VIEWER_STDOUT_WORLD);

    // int nexten;
    // MatGetLocalSize(JJJ, &nexten, &nexten);

    // int isna_arr[4] = {8, 4, 6, 0};
    // IS isna, isna1;
    // // ISCreateStride(PETSC_COMM_SELF, 1, 4, 1, &isna);
    // ISCreateGeneral(PETSC_COMM_WORLD, 1 + 1*(rank==2), isna_arr+rank, PETSC_COPY_VALUES, &isna);
    // ISIntersect(isna, isrow1, &isna1);

    // IS isrow2, isrow3, isrow4;
    // ISGlobalToLocalMappingApplyIS(itog, IS_GTOLM_DROP, isrow1, &isrow2);
    // ISLocalToGlobalMappingApplyIS(ltog, isrow2, &isrow3);
    // ISGlobalToLocalMappingApplyIS(ltog, IS_GTOLM_DROP, isna, &isrow4);

    // ISView(isna, PETSC_VIEWER_STDOUT_WORLD);
    // ISView(isna1, PETSC_VIEWER_STDOUT_WORLD);
    // // if(rank==0) ISView(isrow2, PETSC_VIEWER_STDOUT_SELF);
    // // if(rank==0) ISView(isrow3, PETSC_VIEWER_STDOUT_SELF);
    // if(rank==0) ISView(isrow4, PETSC_VIEWER_STDOUT_SELF);

    // int islistlen;
    // IS* islist;
    // ISPairToList(isrow, isrow1, &islistlen, islist);
    // // ISListToPair(PETSC_COMM_WORLD, 4, &isna, &islist1, islist2);
    // // ISView(islist1, PETSC_VIEWER_STDOUT_WORLD);
    // // ISView(islist2, PETSC_VIEWER_STDOUT_WORLD);
    // ISView(islist[0], PETSC_VIEWER_STDOUT_WORLD);


    // // MatGetOwnershipRange(J2, &istart, &iend);
    // // printf("Rank %i:\t%i %i\n", rank, istart, iend);
    // // MatGetOwnershipRangeColumn(J2, &istart, &iend);
    // // printf("Rank %i:\t%i %i\n", rank, istart, iend);

    // // VecStrideSubSetGather


    // // /* Set dimensions */
    // // if(verbose) printf("\tSetting dimensions...\n");
    // // int nb;
    // // int nt, ns, nu;
    // // int mt, ms, mu;
    // // MatGetSizeLoad("data/At", &mt, &nt);
    // // MatGetSizeLoad("data/As", &ms, &ns);
    // // MatGetSizeLoad("data/Ab", &mu, &nb);
    // // nu = nt * ns;
    // // mu = mt * ms;
        
    // // int nt_local, nt_ghost, nt_exten, nt_separ;
    // // nt_local = graph->offset[1] - graph->offset[0];
    // // nt_ghost = graph->offset[2] - graph->offset[1];
    // // nt_exten = graph->offset[2] - graph->offset[0];
    // // nt_separ = graph->offset[3] - graph->offset[2];
    
    // // int nu_local, nu_ghost, nu_exten, nu_separ;
    // // nu_local = nt_local * ns;
    // // nu_ghost = nt_ghost * ns;
    // // nu_exten = nt_exten * ns;
    // // nu_separ = nt_separ * ns;

    // // int mt_local, mt_ghost, mt_exten, mt_separ;
    // // mt_local = nt_local;
    // // mt_ghost = nt_ghost;
    // // mt_exten = nt_exten;
    // // mt_separ = nt_separ;

    // // int mu_local, mu_ghost, mu_exten, mu_separ;
    // // mu_local = mt_local * ms;
    // // mu_ghost = mt_ghost * ms;
    // // mu_exten = mt_exten * ms;
    // // mu_separ = mt_separ * ms;

    // // if(verbose) printf("\tObservation size:\tms=%i, mt=%i, mu=%i, nb=%i\n", ms, mt, mu, nb);
    // // if(verbose) printf("\tLatent field size:\tns=%i, nt=%i, nu=%i\n", ns, nt, nu);
    // // if(verbose) printf("\tLocal temporal size:\tnt_local=%i, nt_ghost=%i, nt_exten=%i, nt_separ=%i\n", nt_local, nt_ghost, nt_exten, nt_separ);
    // // if(verbose) printf("\tLocal latent size:\tnu_local=%i, nu_ghost=%i, nu_exten=%i, nu_separ=%i\n", nu_local, nu_ghost, nu_exten, nu_separ);


    // // /* Load the data */
    // // if(verbose) printf("\tLoading the data...\n");
    // // Mat J0, J1, J2, At;
    // // Mat K1, K2, K3, As;
    // // Mat J0_seq, J1_seq, J2_seq, At_seq;
    // // Mat Ab;
    // // Vec y;
    // // IS is_na;

    // // MatCreateLoad(PETSC_COMM_WORLD, MATMPIAIJ, nt_local, nt_local, PETSC_DETERMINE, PETSC_DETERMINE, "data/J0", &J0);
    // // MatCreateLoad(PETSC_COMM_WORLD, MATMPIAIJ, nt_local, nt_local, PETSC_DETERMINE, PETSC_DETERMINE, "data/J1", &J1);
    // // MatCreateLoad(PETSC_COMM_WORLD, MATMPIAIJ, nt_local, nt_local, PETSC_DETERMINE, PETSC_DETERMINE, "data/J2", &J2);
    // // MatCreateLoad(PETSC_COMM_WORLD, MATMPIAIJ, nt_local, nt_local, PETSC_DETERMINE, PETSC_DETERMINE, "data/At", &At);
    // // MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/K1", &K1);
    // // MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/K2", &K2);
    // // MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/K3", &K3);
    // // MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/As", &As);
    // // MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/J0", &J0_seq);
    // // MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/J1", &J1_seq);
    // // MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/J2", &J2_seq);
    // // MatCreateLoad(PETSC_COMM_SELF, MATSEQAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/At", &At_seq);
    // // MatCreateLoad(PETSC_COMM_WORLD, MATMPIDENSE, mu_local, nb*last, PETSC_DETERMINE, nb, "data/Ab", &Ab);
    // // VecCreateLoad(PETSC_COMM_WORLD, VECMPI, mu_local, PETSC_DETERMINE, "data/y", &y);
    // // ISCreateLoad(PETSC_COMM_SELF, "data/isna", &is_na);


    // // /* Profiling checkpoint */
    // // if(profile) PetscTime(&t_end);
    // // if(profile) PetscMallocGetCurrentUsage(&mem_end);
    // // if(verbose) printf("\tTime spent:\t\t%f sec\n", t_end - t_start);
    // // if(verbose) printf("\tMemory allocated:\t%f bytes\n", mem_end - mem_start);


    // // /* ---------------------------------------------------------------- */
    // // /* ---------------------------- ASSEMBLY -------------------------- */
    // // /* ---------------------------------------------------------------- */
    // // if(verbose) printf("\nASSEMBLY\n");
    // // if(profile) PetscTime(&t_start);
    // // if(profile) PetscMallocGetCurrentUsage(&mem_start);


    // // /* Create index sets */
    // // if(verbose) printf("\tCreating index sets...\n");
    // // IS isnt_local, isnt_ghost, isnt_exten, isnt_separ;
    // // IS isnu_local, isnu_ghost, isnu_exten, isnu_separ;
    // // IS ismt_local, ismt_ghost, ismt_exten, ismt_separ;
    // // IS ismu_local, ismu_ghost, ismu_exten, ismu_separ;
    // // ISCreateGeneral(PETSC_COMM_SELF, nt_local, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &isnt_local);
    // // ISCreateGeneral(PETSC_COMM_SELF, nt_ghost, graph->vert + graph->offset[1], PETSC_COPY_VALUES, &isnt_ghost);
    // // ISCreateGeneral(PETSC_COMM_SELF, nt_exten, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &isnt_exten);
    // // ISCreateGeneral(PETSC_COMM_SELF, nt_separ, graph->vert + graph->offset[2], PETSC_COPY_VALUES, &isnt_separ);
    // // ISCreateBlock(PETSC_COMM_SELF, ns, nt_local, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &isnu_local);
    // // ISCreateBlock(PETSC_COMM_SELF, ns, nt_ghost, graph->vert + graph->offset[1], PETSC_COPY_VALUES, &isnu_ghost);
    // // ISCreateBlock(PETSC_COMM_SELF, ns, nt_exten, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &isnu_exten);
    // // ISCreateBlock(PETSC_COMM_SELF, ns, nt_separ, graph->vert + graph->offset[2], PETSC_COPY_VALUES, &isnu_separ);
    // // ISCreateGeneral(PETSC_COMM_SELF, mt_local, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &ismt_local);
    // // ISCreateGeneral(PETSC_COMM_SELF, mt_ghost, graph->vert + graph->offset[1], PETSC_COPY_VALUES, &ismt_ghost);
    // // ISCreateGeneral(PETSC_COMM_SELF, mt_exten, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &ismt_exten);
    // // ISCreateGeneral(PETSC_COMM_SELF, mt_separ, graph->vert + graph->offset[2], PETSC_COPY_VALUES, &ismt_separ);
    // // ISCreateBlock(PETSC_COMM_SELF, ms, mt_local, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &ismu_local);
    // // ISCreateBlock(PETSC_COMM_SELF, ms, mt_ghost, graph->vert + graph->offset[1], PETSC_COPY_VALUES, &ismu_ghost);
    // // ISCreateBlock(PETSC_COMM_SELF, ms, mt_exten, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &ismu_exten);
    // // ISCreateBlock(PETSC_COMM_SELF, ms, mt_separ, graph->vert + graph->offset[2], PETSC_COPY_VALUES, &ismu_separ);
    // // // ISCreateStride(PETSC_COMM_SELF, nu_local, 0, 1, &is_local2);
    // // // Check out VecISAXPY...


    // // /* Kronecker products */
    // // if(verbose) printf("\tAssembling the precision matrix...\n");
    // // Mat Quu, Quu1, Quu2, Quu3, Au;
    // // Mat Qub;
    // // Mat Qbb, Qbb1;
    // // MatMPIAIJKron(J0, K3, &Quu);
    // // MatMPIAIJKron(J1, K2, &Quu1);
    // // MatMPIAIJKron(J2, K1, &Quu2);
    // // MatMPIAIJKron(At, As, &Au);
    // // MatZeroRowsIS(Ab, is_na, 0.0, y, y);
    // // MatZeroRowsIS(Au, is_na, 0.0, NULL, NULL);
    // // MatCreateDense(PETSC_COMM_WORLD, nb*last, nb*last, nb, nb, NULL, &Qbb);

    // // MatTransposeMatMult(Au, Au, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Quu3);
    // // MatAXPY(Quu, 1.0, Quu1, DIFFERENT_NONZERO_PATTERN);
    // // MatAXPY(Quu, 1.0, Quu2, DIFFERENT_NONZERO_PATTERN);
    // // MatAXPY(Quu, 1.0, Quu3, DIFFERENT_NONZERO_PATTERN);
    // // MatTransposeMatMult(Ab, Ab, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Qbb1);
    // // MatAXPY(Qbb, 1.0, Qbb1, SAME_NONZERO_PATTERN);
    // // MatTransposeMatMult(Au, Ab, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Qub);


    // // /* Assemble the extended local matrix */
    // // if(verbose) printf("\tAssembling the extended local precision...\n");
    // // Mat J0_exten, J1_exten, J2_exten, At_exten;
    // // Mat Quu_exten, Quu_exten1, Quu_exten2, Au_exten;
    // // MatCreateSubMatrix(J0_seq, isnt_exten, isnt_exten, MAT_INITIAL_MATRIX, &J0_exten);
    // // MatCreateSubMatrix(J1_seq, isnt_exten, isnt_exten, MAT_INITIAL_MATRIX, &J1_exten);
    // // MatCreateSubMatrix(J2_seq, isnt_exten, isnt_exten, MAT_INITIAL_MATRIX, &J2_exten);
    // // MatCreateSubMatrix(At_seq, ismt_exten, isnt_exten, MAT_INITIAL_MATRIX, &At_exten);
    // // MatSeqAIJKron(J0_exten, K3, MAT_INITIAL_MATRIX, &Quu_exten);
    // // MatSeqAIJKron(J1_exten, K2, MAT_INITIAL_MATRIX, &Quu_exten1);
    // // MatSeqAIJKron(J2_exten, K1, MAT_INITIAL_MATRIX, &Quu_exten2);
    // // MatSeqAIJKron(At_exten, As, MAT_INITIAL_MATRIX, &Au_exten);
    

    // // // if(!rank) ISView(is_na, PETSC_VIEWER_STDOUT_SELF);
    // // // IS is_na_exten;
    // // // ISIntersect(is_na, ismu_exten, &is_na_exten);
    // // // if(!rank) ISView(ismu_exten, PETSC_VIEWER_STDOUT_SELF);
    // // // if(!rank) ISView(is_na_exten, PETSC_VIEWER_STDOUT_SELF);
    // // // MatZeroRowsIS(Au_exten, is_na_exten, 0.0, NULL, NULL);
    // // // if(!rank) MatView(Au_exten, PETSC_VIEWER_STDOUT_SELF);


    // // /* !!!!!!!!!!!!!!!         LOCAL/GLOBAL MAPPING        !!!!!!!!!!!!!!!!!! */
    // // /* !!!!!!!!!!!!!!!         LOCAL/GLOBAL MAPPING        !!!!!!!!!!!!!!!!!! */
    // // /* !!!!!!!!!!!!!!!         LOCAL/GLOBAL MAPPING        !!!!!!!!!!!!!!!!!! */
    // // // IS is_na_exten;
    // // // ISLocalToGlobalMapping ltog_na_exten;
    // // // ISLocalToGlobalMappingCreateIS(ismu_exten, &ltog_na_exten);
    // // // ISLocalToGlobalMappingSetType(ltog_na_exten, ISLOCALTOGLOBALMAPPINGHASH);
    // // // ISGlobalToLocalMappingApplyIS(ltog_na_exten, IS_GTOLM_DROP, is_na, &is_na_exten);
    // // // MatZeroRowsIS(Au_exten, is_na_exten, 0.0, NULL, NULL);
    // // // if(!rank) ISView(ismu_local, PETSC_VIEWER_STDOUT_SELF);
    // // // if(!rank) ISView(is_na, PETSC_VIEWER_STDOUT_SELF);
    // // // if(!rank) ISView(is_na_exten, PETSC_VIEWER_STDOUT_SELF);
    // // // if(!rank) MatView(Au_exten, PETSC_VIEWER_STDOUT_SELF);
    // // // if(!rank) ISLocalToGlobalMappingView(ltog_na_exten, PETSC_VIEWER_STDOUT_SELF);


    // // /* !!!!!!!!!!!!!!!         OVERLAPPING        !!!!!!!!!!!!!!!!!! */
    // // /* !!!!!!!!!!!!!!!         OVERLAPPING        !!!!!!!!!!!!!!!!!! */
    // // /* !!!!!!!!!!!!!!!         OVERLAPPING        !!!!!!!!!!!!!!!!!! */
    // // // Mat *JJJ[1];
    // // // IS isnt_isnt;
    // // // ISCreate(PETSC_COMM_WORLD, &isnt_isnt);
    // // // ISCreateGeneral(PETSC_COMM_WORLD, nt_local, graph->vert + graph->offset[0], PETSC_COPY_VALUES, &isnt_isnt);
    // // // // ISCreateStride(PETSC_COMM_WORLD, nt, 0, 1, &isnt_isnt);
    // // // // MatGetOwnershipIS(J2, &isnt_isnt, NULL);
    // // // IS isis[1];
    // // // isis[0] = isnt_isnt;
    // // // MatView(J2, PETSC_VIEWER_STDOUT_WORLD);
    // // // ISView(isis[0], PETSC_VIEWER_STDOUT_WORLD);
    // // // MatCreateSubMatrices(J2, 1, isis, isis, MAT_INITIAL_MATRIX, JJJ);
    // // // if(!rank) MatView(*JJJ[0], PETSC_VIEWER_STDOUT_SELF);
    // // // MatIncreaseOverlap(J2, 1, isis, 1);     //  does not ignore zeros
    // // // MatCreateSubMatrices(J2, 1, isis, isis, MAT_INITIAL_MATRIX, JJJ);
    // // // ISView(isis[0], PETSC_VIEWER_STDOUT_WORLD);
    // // // if(!rank) MatView(*JJJ[0], PETSC_VIEWER_STDOUT_SELF);


    // // /* !!!!!!!!!!!!!!!         PARTITIONING        !!!!!!!!!!!!!!!!!! */
    // // /* !!!!!!!!!!!!!!!         PARTITIONING        !!!!!!!!!!!!!!!!!! */
    // // /* !!!!!!!!!!!!!!!         PARTITIONING        !!!!!!!!!!!!!!!!!! */
    // // Mat KK, KKK, KKKK;
    // // MatCreateLoad(PETSC_COMM_WORLD, MATMPIAIJ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE, "data/J2", &KK);
    // // int part_count[size];
    // // IS is_part, is_num, is_perm, is_permi;
    // // MatPartitioning part;
    // // MatPartitioningCreate(PETSC_COMM_WORLD, &part);
    // // MatPartitioningSetAdjacency(part, KK);
    // // MatPartitioningSetNParts(part, size);
    // // MatPartitioningSetType(part, MATPARTITIONINGPARMETIS);
    // // MatPartitioningApply(part, &is_part);
    // // ISPartitioningToNumbering(is_part, &is_num);
    // // ISPartitioningCount(is_part, size, part_count);
    // // ISInvertPermutation(is_num, part_count[rank], &is_perm);
    // // MatCreateSubMatrix(KK, is_perm, is_perm, MAT_INITIAL_MATRIX, &KKK);
    // // // MatPermute(KK, is_perm, is_perm, &KKK);
    // // ISInvertPermutation(is_perm, PETSC_DECIDE, &is_permi);
    // // MatCreateSubMatrix(KKK, is_permi, is_permi, MAT_INITIAL_MATRIX, &KKKK);
    // // // MatPermute(KKK, is_permi, is_permi, &KKKK);
    

    // // ISView(is_part, PETSC_VIEWER_STDOUT_WORLD);
    // // ISView(is_num, PETSC_VIEWER_STDOUT_WORLD);
    // // ISView(is_perm, PETSC_VIEWER_STDOUT_WORLD);
    // // ISView(is_permi, PETSC_VIEWER_STDOUT_WORLD);
    // // // MatView(KK, PETSC_VIEWER_STDOUT_WORLD);
    // // // MatView(KKKK, PETSC_VIEWER_STDOUT_WORLD);

    // // IS is_over;
    // // ISDuplicate(is_permi, &is_over);
    // // ISCopy(is_permi, is_over);
    // // MatIncreaseOverlap(KKKK, 1, &is_over, 1);
    // // ISView(is_over, PETSC_VIEWER_STDOUT_WORLD);

    // // IS is_best;
    // // int part_cumsum[size];
    // // part_cumsum[0] = 0;
    // // for(int i=1; i<size; i++) part_cumsum[i] = part_count[i-1] + part_cumsum[i-1];
    // // if(!rank) for(int i=0; i<size; i++) printf("%i %i\n", part_count[i], part_cumsum[i]);
    // // ISCreateStride(PETSC_COMM_WORLD, part_count[rank], part_cumsum[rank], 1, &is_best);
    // // ISToGeneral(is_best);
    // // ISView(is_best, PETSC_VIEWER_STDOUT_WORLD);

    // // Mat SS;
    // // MatCreateSubMatrixSeq(KK, is_best, is_best, &SS);
    // // MatView(KK, PETSC_VIEWER_STDOUT_WORLD);
    // // if(rank==2) MatView(SS, PETSC_VIEWER_STDOUT_SELF);
    // // MatCreateSubMatrix(KK, is_best, is_best, MAT_INITIAL_MATRIX, &SS);

    // // int istart, iend;
    // // MatGetOwnershipRange(KK, &istart, &iend);
    // // printf("Rank %i:\t%i %i\n", rank, istart, iend);
    // // MatGetOwnershipRange(SS, &istart, &iend);
    // // printf("Rank %i:\t%i %i\n", rank, istart, iend);