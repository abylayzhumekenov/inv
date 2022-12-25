#include "inv_matrix.h"


/**
 * @brief Compute a Kronecker product of distributed and sequential matrices
 * 
 * @param A An MPI matrix
 * @param B A Seq matrix
 * @param reuse Use MAT_INTIAL_MATRIX for cleaner memory management
 * @param C An MPI Kronecker product
 * @return PetscErrorCode 
 */
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