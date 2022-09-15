#include "inv_solver.h"


KSP solver_ksp(MPI_Comm comm, Mat Q, int max_niter, int verbose){

    /* Create a GMRES solver */
    if(verbose) printf("%sSetting up a solver..\n", INV_SOLVER_VERBOSE);
    KSP ksp;
    PC pc;
    KSPCreate(comm, &ksp);
    KSPSetOperators(ksp, Q, Q);
    KSPSetType(ksp, KSPGMRES);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCICC);  // Later fix to some parallel preconditioner
    PCSetFromOptions(pc);
    // KSPSetPCSide(ksp, PC_SYMMETRIC);
    // KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, max_niter);
    // KSPSetComputeEigenvalues(ksp, PETSC_TRUE);
    // KSPSetComputeRitz(ksp, PETSC_TRUE);
    // KSPGMRESSetRestart(ksp, max_niter);
    // KSPGMRESSetOrthogonalization(ksp,KSPGMRESClassicalGramSchmidtOrthogonalization);
    // KSPGMRESSetCGSRefinementType(ksp, KSP_GMRES_CGS_REFINE_IFNEEDED);                   // Without this, even CGS fails...
    // KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization);  // MGS fails too...
    KSPSetFromOptions(ksp);
    if(verbose) printf("%sSolver created.\n", INV_SOLVER_VERBOSE);

    return ksp;
}


Vec solver_sample_contribution(KSP ksp, Mat Q, Vec x, IS is_x, IS is_a, int verbose){

    if(verbose) printf("%sComputing the sample contribution..\n", INV_SOLVER_VERBOSE);
    Vec y, w;
    
    /* Compute Q[a,s] %*% x[s] */
    VecISSet(x, is_x, 0.0);
    VecDuplicate(x, &y);
    MatMult(Q, x, y);
    VecGetSubVector(y, is_a, &y);
    VecDuplicate(y, &w);

    /* Solve with Q[a,a] and square */
    KSPSolve(ksp, y, w);
    VecPointwiseMult(w, w, w);
    if(verbose) printf("%sSample contribution computed.\n", INV_SOLVER_VERBOSE);

    return w;
}


Vec solver_base_case(MPI_Comm comm, Mat Q, int verbose){

    int n;
    Vec d;
    Mat L, II, S;
    IS isrow, iscol;
    
    MatGetSize(Q, &n, &n);
    VecCreate(PETSC_COMM_SELF, &d);
    VecSetType(d, VECSEQ);
    VecSetSizes(d, PETSC_DECIDE, n);
    VecSetFromOptions(d);
    VecSetUp(d);

    if(verbose) printf("%sSolving the base case.\n", INV_SOLVER_VERBOSE);
    MatSetOption(Q, MAT_SYMMETRIC, PETSC_TRUE);
    MatGetOrdering(Q, MATORDERINGNATURAL, &isrow, &iscol);
    MatDuplicate(Q, MAT_SHARE_NONZERO_PATTERN, &L);
    MatGetFactor(Q, MATSOLVERMUMPS, MAT_FACTOR_CHOLESKY, &L);
    MatCholeskyFactorSymbolic(L, Q, isrow, NULL);
    MatCholeskyFactorNumeric(L, Q, NULL);
    MatCreateDense(comm, PETSC_DECIDE, PETSC_DECIDE, n, n, NULL, &II);
    MatCreateDense(comm, PETSC_DECIDE, PETSC_DECIDE, n, n, NULL, &S);
    for(int i=0; i<n; i++){
        MatSetValue(II, i, i, 1.0, INSERT_VALUES);
    }
    MatAssemblyBegin(II, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(II, MAT_FINAL_ASSEMBLY);
    MatMatSolve(L, II, S);
    MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
    MatGetDiagonal(S, d);
    if(verbose) printf("%sBase case solved.\n", INV_SOLVER_VERBOSE);

    MatDestroy(&L);
    MatDestroy(&II);
    MatDestroy(&S);
    ISDestroy(&isrow);
    ISDestroy(&iscol);
    if(verbose) printf("%sClean up successful.\n", INV_SOLVER_VERBOSE);

    return d;
}


Vec solver_assemble_solution(Vec d, IS is_part, IS is_a, int n, int verbose){

    if(verbose) printf("%sAssembling local solution..\n", INV_SOLVER_VERBOSE);
    Vec vv, v;
    VecCreate(PETSC_COMM_SELF, &vv);
    VecSetType(vv, VECSEQ);
    VecSetSizes(vv, PETSC_DECIDE, n);
    VecSetFromOptions(vv);
    VecSetUp(vv);
    VecSet(vv, 0.0);
    VecISCopy(vv, is_a, SCATTER_FORWARD, d);
    VecGetSubVector(vv, is_part, &vv);

    if(verbose) printf("%sAssembling global solution..\n", INV_SOLVER_VERBOSE);
    VecCreate(PETSC_COMM_WORLD, &v);
    VecSetType(v, VECMPI);
    VecSetSizes(v, PETSC_DECIDE, n);
    VecSetFromOptions(v);
    VecSetUp(v);

    int n_is;
    double* vv_temp;
    const int* is_part_temp;
    VecGetArray(vv, &vv_temp);
    ISGetSize(is_part, &n_is);
    ISGetIndices(is_part, &is_part_temp);
    for(int i=0; i<n_is; i++){
        VecSetValues(v, n_is, is_part_temp, vv_temp, INSERT_VALUES);
    }
    VecAssemblyBegin(v);
    VecAssemblyEnd(v);
    VecRestoreArray(vv, &vv_temp);
    VecDestroy(&vv);
    vv_temp = NULL;
    is_part_temp = NULL;
    if(verbose) printf("%sAssembling complete.\n", INV_SOLVER_VERBOSE);

    return v;
}