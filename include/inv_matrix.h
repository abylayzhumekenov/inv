#ifndef INV_MATRIX_H
#define INV_MATRIX_H

#include <petsc.h>

#include "inv_graph.h"
#include "inv_is.h"


Mat InvMatrixCreateFromFile(MPI_Comm comm, char* filename, InvIS* mapping, int verbose);


Mat InvMatrixCreateLocalFromFile(MPI_Comm comm, char* filename, InvIS* mapping, int verbose);


Mat InvMatrixFactorLocal(Mat Q, int verbose);


Mat InvMatrixCreateLocalSubmatrix(MPI_Comm comm, Mat Q, InvIS* mapping, InvIS* mapping_sub, int verbose);


Mat InvMatrixCreateRW2(MPI_Comm comm, int n, double alpha, double tau, int cyclic, int verbose);


#endif