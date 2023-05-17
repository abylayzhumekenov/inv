#include <petsc.h>


/**
 * @brief Load a matrix from a file
 * 
 * @param comm 
 * @param type 
 * @param m 
 * @param n 
 * @param M 
 * @param N 
 * @param filename 
 * @param A 
 * @return PetscErrorCode 
 */
PetscErrorCode MatCreateLoad(MPI_Comm comm, MatType type, PetscInt m, PetscInt n, PetscInt M, PetscInt N, const char filename[], Mat* A);


/**
 * @brief Return a matrix size from a file header without loading
 * 
 * @param filename 
 * @param m 
 * @param n 
 * @return PetscErrorCode 
 */
PetscErrorCode MatGetSizeLoad(const char filename[], int* m, int* n);


/**
 * @brief Create a sequential submatrix defined by local portions of index sets
 * 
 * @param A 
 * @param isrow 
 * @param iscol 
 * @param B 
 * @return PetscErrorCode 
 */
PetscErrorCode MatCreateSubMatrixSeq(Mat A, IS isrow, IS iscol, Mat* B);


/**
 * @brief Compute the Kronecker product of a parallel and a sequential matrices
 * 
 * @param A 
 * @param B 
 * @param C 
 * @return PetscErrorCode 
 */
PetscErrorCode MatMPIAIJKron(Mat A, Mat B, Mat* C);


/**
 * @brief Invert a dense matrix in place
 * 
 * @param A 
 * @return PetscErrorCode 
 */
PetscErrorCode MatDenseInvertLapack(Mat A);


/**
 * @brief Compute an element wise product of two dense matrices in place
 * 
 * @param A 
 * @param B 
 * @return PetscErrorCode 
 */
PetscErrorCode MatDensePointwiseMult(Mat A, Mat B);


/**
 * @brief Compute selected elements of the inverse using MUMPS
 * 
 * @param A 
 * @param B 
 * @return PetscErrorCode 
 */
PetscErrorCode MatSeqAIJInvertMUMPS(Mat A, Mat* B, int verbose);


/**
 * @brief Compute selected elements of the inverse using PARDISO
 * 
 * @param A 
 * @param B 
 * @return PetscErrorCode 
 */
PetscErrorCode MatSeqAIJInvertPARDISO(Mat A, Mat* B, int verbose);