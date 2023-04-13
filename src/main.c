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
    MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- PARTITIONING ---------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nPARTITIONING\n");
    MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- ASSEMBLY -------------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nASSEMBLY\n");
    MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));

    
    /* ---------------------------------------------------------------- */
    /* ---------------------------- DIRECT SOLVE ---------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nDIRECT SOLVE\n");
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Factor the extended matrix */
    if(verbose) printf("\tFactoring the extended matrix...\n");
    Mat L_exten;
    IS is_factor;
    MatSetOption(Quu_exten, MAT_SYMMETRIC, PETSC_TRUE);
    MatGetOrdering(Quu_exten, MATORDERINGNATURAL, &is_factor, &is_factor);
    MatGetFactor(Quu_exten, MATSOLVERMUMPS, MAT_FACTOR_CHOLESKY, &L_exten);
    MatCholeskyFactorSymbolic(L_exten, Quu_exten, is_factor, NULL);
    MatCholeskyFactorNumeric(L_exten, Quu_exten, NULL);


    /* Selected inversion */
    if(verbose) printf("\tSelected inversion...\n");
    Vec d_exten;
    Mat D_exten;
    VecCreateSeq(PETSC_COMM_SELF, nu_exten, &d_exten);
    MatCreate(PETSC_COMM_SELF, &D_exten);
    MatSetSizes(D_exten, nu_exten, nu_exten, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetType(D_exten, MATSEQAIJ);
    MatSetFromOptions(D_exten);
    MatSetUp(D_exten);
    for(int i=0; i<nu_exten; i++) MatSetValue(D_exten, i, i, 1.0, INSERT_VALUES);
    MatAssemblyBegin(D_exten, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(D_exten, MAT_FINAL_ASSEMBLY);
    MatMumpsGetInverseTranspose(L_exten, D_exten);
    MatGetDiagonal(D_exten, d_exten);


    /* Profiling checkpoint */
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


        // /* ---------------------------------------------------------------- */
        // /* ---------------------------- CHECK PURELY DIRECT SOLVE --------- */
        // /* ---------------------------------------------------------------- */
        // if(verbose) printf("\nCHECK PURELY DIRECT SOLVE\n");
        // if(profile) PetscTime(&t_start);
        // if(profile) PetscMallocGetCurrentUsage(&mem_start);

        // Mat LL;
        // IS isis;
        // MatSetOption(Quu, MAT_SYMMETRIC, PETSC_TRUE);
        // MatGetOrdering(Quu, MATORDERINGNATURAL, &isis, &isis);
        // MatGetFactor(Quu, MATSOLVERMUMPS, MAT_FACTOR_CHOLESKY, &LL);
        // MatCholeskyFactorSymbolic(LL, Quu, isis, NULL);
        // MatCholeskyFactorNumeric(LL, Quu, NULL);

        // /* Selected inversion */
        // if(verbose) printf("\tSelected inversion...\n");
        // Vec dd;
        // Mat DD;
        // int istart, iend;
        // VecCreateMPI(PETSC_COMM_WORLD, nu*(!rank), PETSC_DETERMINE, &dd);
        // MatCreate(PETSC_COMM_WORLD, &DD);
        // MatSetSizes(DD, nu*(!rank), nu*(!rank), PETSC_DETERMINE, PETSC_DETERMINE);
        // MatSetType(DD, MATMPIAIJ);
        // MatSetFromOptions(DD);
        // MatSetUp(DD);
        // MatGetOwnershipRange(DD, &istart, &iend);
        // for(int i=istart; i<iend; i++) MatSetValue(DD, i, i, 1.0, INSERT_VALUES);
        // MatAssemblyBegin(DD, MAT_FINAL_ASSEMBLY);
        // MatAssemblyEnd(DD, MAT_FINAL_ASSEMBLY);
        // MatMumpsGetInverseTranspose(LL, DD);
        // MatGetDiagonal(DD, dd);

        // /* Profiling checkpoint */
        // if(profile) PetscTime(&t_end);
        // if(profile) PetscMallocGetCurrentUsage(&mem_end);
        // if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
        // if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- KSP SOLVERS ----------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nKSP SOLVERS\n");
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_start);
    if(profile) PetscMallocGetCurrentUsage(&mem_start);


    /* Create a solver for sampling */
    if(verbose) printf("\tSetting up a sampling solver...\n");
    InvShellPC* shell;
    KSP ksp_sampling;
    PC pc_sampling;
        double tt0, tt1;
        PetscTime(&tt0);
    KSPCreate(PETSC_COMM_WORLD, &ksp_sampling);
        PetscTime(&tt1);
        printf("\n\tTime spent on %i:\t\t%f sec\n", rank, tt1-tt0);
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
        KSPSetPCSide(ksp_sampling, PC_SYMMETRIC);
    KSPGetPC(ksp_correction, &pc_correction);
    PCSetType(pc_correction, PCBJACOBI);
    PCSetFromOptions(pc_correction);
    KSPSetTolerances(ksp_correction, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, n_iter);
    KSPSetFromOptions(ksp_correction);


    /* Create a solver for covariates */
    if(verbose) printf("\tSetting up a covariates solver...\n");
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
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- SAMPLING CORRECTION --------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nSAMPLING CORRECTION\n");
    MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- COVARIATE CORRECTION -------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nCOVARIATE CORRECTION\n");
    MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    /* ---------------------------------------------------------------- */
    /* ---------------------------- SOLVE FOR THE MEAN ---------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nSOLVE FOR THE MEAN\n");
    MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));

        int nnn;
        const char* reason;
        double sumy;
        VecSum(y, &sumy);
        KSPGetIterationNumber(ksp_covariates, &nnn);
        KSPGetConvergedReasonString(ksp_covariates, &reason);
        if(profile) printf("\t\tReason:\t%s\tIterations:\t%i\tSum of y:\t%f\n", reason, nnn, sumy);

        if(profile) PetscTime(&t_start);
        PCApplySymmetricLeft(pc_covariates, bu, meanu);
        if(profile) PetscTime(&t_end);
        if(profile) printf("\n\tUsual PCApplyLeft:\t\t%f sec\n", t_end - t_start);
        if(profile) PetscTime(&t_start);
        PCApplySymmetricRight(pc_covariates, bu, meanu);
        if(profile) PetscTime(&t_end);
        if(profile) printf("\n\tUsual PCApplyRight:\t\t%f sec\n", t_end - t_start);
        if(profile) PetscTime(&t_start);
        MatMult(Quu, bu, meanu);
        if(profile) PetscTime(&t_end);
        if(profile) printf("\n\tUsual MatMult:\t\t%f sec\n", t_end - t_start);


    /* ---------------------------------------------------------------- */
    /* ---------------------------- CLEAN UP -------------------------- */
    /* ---------------------------------------------------------------- */
    if(verbose) printf("\nCLEAN UP\n");
    MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Barrier(MPI_COMM_WORLD);
    if(profile) PetscTime(&t_end);
    if(profile) PetscMallocGetCurrentUsage(&mem_end);
    if(profile) printf("\n\tTime spent:\t\t%f sec\n", t_end - t_start);
    if(profile) printf("\tMemory allocated:\t%i bytes\n", (int)(mem_end - mem_start));


    // /* ---------------------------------------------------------------- */
    // /* ---------------------------- FINALIZE -------------------------- */
    // /* ---------------------------------------------------------------- */
    if(verbose) printf("\nFINALIZE\n");


    /* Final profiling */
    MPI_Barrier(MPI_COMM_WORLD);
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
