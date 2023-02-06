#ifndef INV_MATRIX_H
#define INV_MATRIX_H

#include <petsc.h>


/**
 * @brief Compute a Kronecker product of distributed and sequential matrices
 * 
 * @param A An MPI matrix
 * @param B A Seq matrix
 * @param reuse Use MAT_INTIAL_MATRIX for cleaner memory management
 * @param C An MPI Kronecker product
 * @return PetscErrorCode 
 */
PetscErrorCode MatMPIAIJKron(Mat A, Mat B, Mat* C);


/**
 * @brief 
 * 
 * @param Q 
 * @param Q_full 
 * @return PetscErrorCode 
 */
PetscErrorCode MatConcatenateIntercept(Mat Q0, Mat* Q, double tau_y, double tau_b);


#endif