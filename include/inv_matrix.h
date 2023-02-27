#include <petsc.h>


PetscErrorCode MatCreateLoad(MPI_Comm comm, MatType type, PetscInt m, PetscInt n, PetscInt M, PetscInt N, const char filename[], Mat* A);


PetscErrorCode MatGetSizeLoad(const char filename[], int* m, int* n);


PetscErrorCode MatCreateSubMatrixSeq(Mat A, IS isrow, IS iscol, Mat* B);


PetscErrorCode MatMPIAIJKron(Mat A, Mat B, Mat* C);


PetscErrorCode MatDenseInvertLapack(Mat A);


PetscErrorCode MatDensePointwiseMult(Mat A, Mat B);