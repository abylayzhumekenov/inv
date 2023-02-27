#include <petsc.h>


PetscErrorCode VecCreateLoad(MPI_Comm comm, VecType type, PetscInt n, PetscInt N, const char filename[], Vec* v);


PetscErrorCode VecCreateSubVector(Vec x, IS is, Vec* y);


PetscErrorCode VecCreateSubVectorSeq(Vec x, IS is, Vec* y);


PetscErrorCode VecCreateSuperVector(Vec x, Vec* y);