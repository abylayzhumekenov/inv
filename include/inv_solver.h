#ifndef INV_SOLVER_H
#define INV_SOLVER_H

#include <petsc.h>
#include "inv_is.h"


Vec InvSolverBaseCase(Mat Q, int verbose);


// KSP InvSolverCreateKSP(MPI_Comm comm, Mat Q, int max_niter, int verbose);


Vec InvSolverMultQx(MPI_Comm comm, Mat Q, Vec x, InvIS* itol2, InvIS* itol, InvIS* itog, int verbose);


Vec InvSolverSampleSq(Mat L, Vec y, int verbose);


Vec InvSolverAssembleSolution(MPI_Comm comm, Vec d, InvIS* itol, InvIS* itog, InvIS* gtoi, int verbose);


#endif