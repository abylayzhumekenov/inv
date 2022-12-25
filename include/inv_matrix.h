#ifndef INV_MATRIX_H
#define INV_MATRIX_H

#include <petsc.h>


PetscErrorCode MatMPIAIJKron(Mat A, Mat B, MatReuse reuse, Mat* C);


#endif