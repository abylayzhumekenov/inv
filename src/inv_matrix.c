#include <petsc.h>
#include <petscblaslapack.h>

#include "inv_matrix.h"


PetscErrorCode MatCreateLoad(MPI_Comm comm, MatType type, PetscInt m, PetscInt n, PetscInt M, PetscInt N, const char filename[], Mat* A){

    int petsc_decide = (m==PETSC_DECIDE && n==PETSC_DECIDE && M==PETSC_DETERMINE && N==PETSC_DETERMINE);
    PetscViewer viewer;
    MatCreate(comm, A);
    MatSetType(*A, type);
    if(!petsc_decide) MatSetSizes(*A, m, n, M, N);
    PetscViewerBinaryOpen(comm, filename, FILE_MODE_READ, &viewer);
    MatLoad(*A, viewer);

    PetscViewerDestroy(&viewer);
    return 0;
}


PetscErrorCode MatGetSizeLoad(const char filename[], int* m, int* n){

    FILE* file = fopen(filename, "rb");
    unsigned int header[4];
    fread(header, 4, 4, file);
    *m = htobe32(header[1]);
    *n = htobe32(header[2]);
    fclose(file);

    return 0;
}


PetscErrorCode MatCreateSubMatrixSeq(Mat A, IS isrow, IS iscol, Mat* B){

    Mat *BB[1];
    MatCreateSubMatrices(A, 1, &isrow, &iscol, MAT_INITIAL_MATRIX, BB);
    *B = *BB[0];
    
    return 0;
}


PetscErrorCode MatMPIAIJKron(Mat A, Mat B, Mat* C){

    int nt_local, mt_local, ns, ms;
    Mat *A_seq[1], C_seq;
    IS isrow, iscol;

    MatGetOwnershipIS(A, &isrow, &iscol);
    MatCreateSubMatrices(A, 1, &isrow, &iscol, MAT_INITIAL_MATRIX, A_seq);
    MatSeqAIJKron(*A_seq[0], B, MAT_INITIAL_MATRIX, &C_seq);
    MatGetSize(*A_seq[0], &nt_local, &mt_local);
    MatGetSize(B, &ns, &ms);
    MatCreateMPIMatConcatenateSeqMat(PETSC_COMM_WORLD, C_seq, nt_local*ms, MAT_INITIAL_MATRIX, C);

    MatDestroyMatrices(1, A_seq);
    MatDestroy(&C_seq);
    return 0;
}


PetscErrorCode MatDenseInvertLapack(Mat A){

    int rank, size, n;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MatGetLocalSize(A, &n, &n);

    double *a_array, work[n];
    int i_array[n], info;
    MatDenseGetArray(A, &a_array);
    if(rank==size-1){
        dgetrf_(&n, &n, a_array, &n, i_array, &info);
        dgetri_(&n, a_array, &n, i_array, work, &n, &info);
    }
    MatDenseRestoreArray(A, &a_array);

    return 0;
}


PetscErrorCode MatDensePointwiseMult(Mat A, Mat B){

    double *a_array, *b_array;
    int n, m;
    MatGetLocalSize(A, &n, NULL);
    MatGetSize(A, NULL, &m);
    MatDenseGetArray(A, &a_array);
    MatDenseGetArray(B, &b_array);
    for(int i=0; i<n*m; i++) a_array[i] *= b_array[i];
    MatDenseRestoreArray(A, &a_array);
    MatDenseRestoreArray(B, &b_array);

    return 0;
}