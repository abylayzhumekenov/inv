#include "inv_matrix.h"


PetscErrorCode MatMPIAIJKron(Mat A, Mat B, MatReuse reuse, Mat* C){

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    int n_local, n;
    Mat *A_seq[1], C_seq;
    IS isrow[1], iscol[1];
    MatGetOwnershipIS(A, &isrow[0], &iscol[0]);
    MatCreateSubMatrices(A, 1, isrow, iscol, reuse, A_seq);

    MatSeqAIJKron(*A_seq[0], B, reuse, &C_seq);
    MatGetSize(C_seq, &n_local, &n);
    MatCreateMPIMatConcatenateSeqMat(PETSC_COMM_WORLD, C_seq, n_local, reuse, C);

    MatDestroyMatrices(1, A_seq);
    MatDestroy(&C_seq);

    return 0;
}


PetscErrorCode MatConcatenateIntercept(Mat Q0, Mat* Q, double tau){
    
    int rank, size, last;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    last = (rank==size-1);

    int n_local, n;
    MatGetLocalSize(Q0, &n_local, &n_local);
    MatGetSize(Q0, &n, &n);

    Mat M1, M2, M3, Q1, Q2, Q3;
    MatCreateSeqAIJ(PETSC_COMM_SELF, n_local, 1, 1, NULL, &M1);
    MatCreateSeqAIJ(PETSC_COMM_SELF, last, size*n_local, size*n_local, NULL, &M2);
    MatCreateSeqAIJ(PETSC_COMM_SELF, last, 1, 1, NULL, &M3);

    int istart, iend;
    MatGetOwnershipRange(Q0, &istart, &iend);
    for(int i=0; i<n_local; i++) MatSetValue(M1, i, 0, tau, INSERT_VALUES);
    if(last) { for(int i=0; i<n; i++) MatSetValue(M2, 0, i, tau, INSERT_VALUES); }
    if(last) MatSetValue(M3, 0, 0, tau*(n+1), INSERT_VALUES);

    MatCreateMPIMatConcatenateSeqMat(PETSC_COMM_WORLD, M1, last, MAT_INITIAL_MATRIX, &Q1);
    MatCreateMPIMatConcatenateSeqMat(PETSC_COMM_WORLD, M2, n_local, MAT_INITIAL_MATRIX, &Q2);
    MatCreateMPIMatConcatenateSeqMat(PETSC_COMM_WORLD, M3, last, MAT_INITIAL_MATRIX, &Q3);

    Mat QQ[4];
    QQ[0] = Q0;
    QQ[1] = Q1;
    QQ[2] = Q2;
    QQ[3] = Q3;
    MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, QQ, Q);

    MatDestroy(&M1);
    MatDestroy(&M2);
    MatDestroy(&M3);

    return 0;
}