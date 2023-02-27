#include <petsc.h>


PetscErrorCode ISCreateLoad(MPI_Comm comm, const char filename[], IS* is);


PetscErrorCode MatGetOwnershipISRows(MPI_Comm comm, Mat A, IS* is);


PetscErrorCode ISCreateBlockIS(MPI_Comm comm, IS is, int bs, IS* isb);