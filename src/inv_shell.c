#include "inv_shell.h"


PetscErrorCode InvShellCreate(InvShellPC** shell){

    InvShellPC* ctx = (InvShellPC*)malloc(sizeof(InvShellPC));
    ctx->Q = NULL;
    ctx->pc = NULL;
    *shell = ctx;
    
    return 0;
}


PetscErrorCode InvShellSetup(PC pc, Mat Q, MPI_Comm comm){

    InvShellPC* shell;
    PCShellGetContext(pc, &shell);

    KSP *ksp_sub;
    PC pc_core, pc_sub;

    PCCreate(comm, &pc_core);
    PCSetOperators(pc_core, Q, Q);
    PCSetType(pc_core, PCBJACOBI);
    PCSetUp(pc_core);

    PCBJacobiGetSubKSP(pc_core, NULL, NULL, &ksp_sub);
    KSPSetType(*ksp_sub, KSPGMRES);
    KSPGetPC(*ksp_sub, &pc_sub);
    PCSetType(pc_sub, PCICC);
    KSPSetPCSide(*ksp_sub, PC_SYMMETRIC);

    shell->Q = Q;
    shell->pc = pc_core;

    return 0;
}


PetscErrorCode InvShellApply(PC pc, Vec x, Vec y){

        double tt0, tt1;
        PetscTime(&tt0);
    InvShellPC *shell;
    PCShellGetContext(pc, &shell);
    PCApplySymmetricLeft(shell->pc, x, y);
        PetscTime(&tt1);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(!rank) printf("\t\tShellApply:\t%f\n", tt1-tt0);

    return 0;
}


PetscErrorCode InvShellApplyTranspose(PC pc, Vec x, Vec y){

        double tt0, tt1;
        PetscTime(&tt0);
    InvShellPC *shell;
    PCShellGetContext(pc, &shell);
    PCApplySymmetricRight(shell->pc, x, y);
        PetscTime(&tt1);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(!rank) printf("\t\tShellApplyTranspose:\t%f\n", tt1-tt0);

    return 0;
}


PetscErrorCode InvShellApplyBA(PC pc, PCSide side, Vec x, Vec y, Vec work){

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        double tt0, tt1, tt2, tt3, tt4;
        PetscTime(&tt0);
    InvShellPC *shell;
    PCShellGetContext(pc, &shell);

        PetscTime(&tt1);
    PCApplySymmetricRight(shell->pc, x, y);
        PetscTime(&tt2);
    MatMult(shell->Q, y, work);
        PetscTime(&tt3);
    PCApplySymmetricLeft(shell->pc, work, y);
        PetscTime(&tt4);
        tt4 -= tt3;
        tt3 -= tt2;
        tt2 -= tt1;
        tt1 -= tt0;
        if(!rank) printf("\t\tShellApplyBA:\t%f\t%f\t%f\t%f\n", tt1, tt2, tt3, tt4);

    return 0;
}


PetscErrorCode InvShellApplyLeft(PC pc, Vec x, Vec y){

        double tt0, tt1;
        PetscTime(&tt0);
    InvShellPC *shell;
    PCShellGetContext(pc, &shell);
    /*  USE ONLY AS A SAMPLER, CANNOT SOLVE SYSTEMS!!!
        If L^-1 is not applied to z when computing the initial residual, 
        we can construct a subspace from z directly, i.e. K((L^-1*Q*L^-T), z).
        However, if we precondition z as usual, we have to premultiply to 
        get L^-1(Lz) = z, in order to obtain the same subspace. Therefore,
        it is easier to tweak ApplySymmetricLeft function so that it does
        nothing. Note that this function is called only once, and if it is
        called somewhere else, you must go back to original implementation. */
    // PCApplySymmetricLeft(shell->pc_core, x, y);  // when preconditioned rhs is L^-1*L*z = z
    VecCopy(x, y);                               // when preconditioned rhs is (nothing)*z
        PetscTime(&tt1);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(!rank) printf("\t\tShellApplyLeft:\t%f\n", tt1-tt0);

    return 0;
}


PetscErrorCode InvShellApplyRight(PC pc, Vec x, Vec y){

        double tt0, tt1;
        PetscTime(&tt0);
    InvShellPC *shell;
    PCShellGetContext(pc, &shell);
    PCApplySymmetricRight(shell->pc, x, y);
        PetscTime(&tt1);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(!rank) printf("\t\tShellApplyRight:\t%f\n", tt1-tt0);

    return 0;
}


PetscErrorCode InvShellDestroy(PC pc){

    InvShellPC *shell;
    PCShellGetContext(pc, &shell);
    PCDestroy(&shell->pc);
    free(shell);
    
    return 0;
}