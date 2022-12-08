#include "inv_solver.h"

#define INV_SOLVER_VERBOSE "SOLVER:\t\t"


Vec InvSolverBaseCase(Mat L, int verbose){

    int n;
    Vec d;
    Mat II, S;

    MatGetSize(L, &n, &n);
    VecCreate(PETSC_COMM_SELF, &d);
    VecSetType(d, VECSEQ);
    VecSetSizes(d, PETSC_DECIDE, n);
    VecSetFromOptions(d);
    VecSetUp(d);

    if(verbose) printf("%sSolving the base case...\n", INV_SOLVER_VERBOSE);
    MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, n, n, NULL, &II);
    MatCreateDense(PETSC_COMM_SELF, PETSC_DECIDE, PETSC_DECIDE, n, n, NULL, &S);
    for(int i=0; i<n; i++){
        MatSetValue(II, i, i, 1.0, INSERT_VALUES);
    }
    MatAssemblyBegin(II, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(II, MAT_FINAL_ASSEMBLY);
    MatMatSolve(L, II, S);
    MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
    MatGetDiagonal(S, d);
    MatDestroy(&II);
    MatDestroy(&S);
    if(verbose) printf("%sBase case solved.\n", INV_SOLVER_VERBOSE);

    return d;
}


KSP InvSolverCreateKSP(MPI_Comm comm, Mat Q, int max_niter, int verbose){

    if(verbose) printf("%sSetting up a solver...\n", INV_SOLVER_VERBOSE);
    KSP ksp;
    PC pc;
    KSPCreate(comm, &ksp);
    KSPSetOperators(ksp, Q, Q);
    KSPSetType(ksp, KSPGMRES);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCBJACOBI);
    PCSetFromOptions(pc);
    KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, max_niter);
    KSPSetFromOptions(ksp);
    if(verbose) printf("%sSolver created.\n", INV_SOLVER_VERBOSE);

    return ksp;
}


Vec InvSolverMultQx(MPI_Comm comm, Mat Q, Vec x, InvIS* itol2, InvIS* itol, InvIS* itog, int verbose){

    /* Get rank and size */
    int rank, size, error;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* Get local and separator indices */
    if(verbose) printf("%sManaging indices...\n", INV_SOLVER_VERBOSE);
    int n_local = itol2->n_local[rank];
    int n_sub = itol->n_local[rank];
    double zero_sub[n_sub];
    int temp_local[n_local], temp_sub[n_sub], i_local = 0, j_local = 0;
    for(int i=0; i<n_local; i++) temp_local[i] = 0;
    for(int i=0; i<itol2->n_vert; i++){
        if(itol2->part[i]==rank){
            temp_local[i_local] = itog->index[i];
            i_local++;
        }
    }
    i_local = 0;
    for(int i=0; i<n_sub; i++){
        temp_sub[i] = 0;
        zero_sub[i] = 0.0;
    }
    for(int i=0; i<itol2->n_vert; i++){
        if(itol2->part[i]==rank && itol->part[i]==rank){
            temp_sub[i_local] = j_local;
            i_local++;
        }
        j_local += itol2->part[i]==rank;
    }

    /* Share a sample */
    if(verbose) printf("%sCommunicating a sample...\n", INV_SOLVER_VERBOSE);
    int n;
    VecGetSize(x, &n);
    double *sendbuf, recvbuf[n];
    VecGetArray(x, &sendbuf);
    error = MPI_Allgatherv(sendbuf, itog->n_local[rank], MPI_DOUBLE, recvbuf, itog->n_local, itog->offset, MPI_DOUBLE, comm);
    if(error && verbose) printf("%sMPI error code %i.\n", INV_SOLVER_VERBOSE, error);
    VecRestoreArray(x, &sendbuf);

    /* Create a sequential vector and compute Q[a,s]*x[s] */
    if(verbose) printf("%sComputing matrix-vector product...\n", INV_SOLVER_VERBOSE);
    Vec xx, x_local, y, y_local;
    IS is_local, is_sub;
    VecCreateSeqWithArray(PETSC_COMM_SELF, 1, n, recvbuf, &xx);
    ISCreateGeneral(PETSC_COMM_SELF, n_local, temp_local, PETSC_COPY_VALUES, &is_local);
    ISCreateGeneral(PETSC_COMM_SELF, n_sub, temp_sub, PETSC_COPY_VALUES, &is_sub);
    VecGetSubVector(xx, is_local, &x_local);
    VecSetValues(x_local, n_sub, temp_sub, zero_sub, INSERT_VALUES);
    VecDuplicate(x_local, &y_local);
    MatMult(Q, x_local, y_local);
    VecGetSubVector(y_local, is_sub, &y);

    VecDestroy(&xx);
    VecDestroy(&x_local);
    VecDestroy(&y_local);
    ISDestroy(&is_local);
    ISDestroy(&is_sub);

    return y;
}


Vec InvSolverSampleSq(Mat L, Vec y, int verbose){

    if(verbose) printf("%sSolving a system...\n", INV_SOLVER_VERBOSE);
    Vec w;
    VecDuplicate(y, &w);
    MatSolve(L, y, w);
    VecPointwiseMult(w, w, w);

    return w;
}


Vec InvSolverAssembleSolution(MPI_Comm comm, Vec d, InvIS* itol, InvIS* itog, InvIS* gtoi, int verbose){

    /* Get rank and size */
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* Compute indices */
    if(verbose) printf("%sManaging indices...\n", INV_SOLVER_VERBOSE);
    int temp[itog->n_local[rank]], i_local = 0, j_local = 0;
    for(int i=0; i<itog->n_vert; i++){
        if(itog->part[i]==rank && itol->part[i]==rank){
            temp[i_local] = j_local;
            i_local++;
        }
        j_local += itol->part[i]==rank;
    }

    /* Assemble the solution */
    if(verbose) printf("%sAssembling the solution...\n", INV_SOLVER_VERBOSE);
    Vec u, v;
    IS is;
    double *u_array;
    VecCreate(comm, &v);
    VecSetType(v, VECMPI);
    VecSetSizes(v, itog->n_local[rank], PETSC_DETERMINE);
    VecSetFromOptions(v);
    VecSetUp(v);

    ISCreateGeneral(PETSC_COMM_SELF, itog->n_local[rank], temp, PETSC_COPY_VALUES, &is);
    VecGetSubVector(d, is, &u);
    VecGetArray(u, &u_array);
    for(int i=0; i<itog->n_local[rank]; i++){
        VecSetValue(v, gtoi->index[i+itog->offset[rank]], u_array[i], INSERT_VALUES);
    }
    VecAssemblyBegin(v);
    VecAssemblyEnd(v);
    VecRestoreArray(u, &u_array);
    VecDestroy(&u);
    if(verbose) printf("%sAssembling complete.\n", INV_SOLVER_VERBOSE);

    return v;
}