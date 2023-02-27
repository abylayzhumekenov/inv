#include <petsc.h>


/**
 * @brief Create a vector from a file
 * 
 * @param comm 
 * @param type 
 * @param n 
 * @param N 
 * @param filename 
 * @param v 
 * @return PetscErrorCode 
 */
PetscErrorCode VecCreateLoad(MPI_Comm comm, VecType type, PetscInt n, PetscInt N, const char filename[], Vec* v);


/**
 * @brief Create a subvector defined by an index set
 * 
 * @param x 
 * @param is 
 * @param y 
 * @return PetscErrorCode 
 */
PetscErrorCode VecCreateSubVector(Vec x, IS is, Vec* y);


/**
 * @brief Create a sequential subvector defined by a local porion of an index set
 * 
 * @param x 
 * @param is 
 * @param y 
 * @return PetscErrorCode 
 */
PetscErrorCode VecCreateSubVectorSeq(Vec x, IS is, Vec* y);


/**
 * @brief Create a parallel vector from sequential ones
 * 
 * @param x 
 * @param y 
 * @return PetscErrorCode 
 */
PetscErrorCode VecCreateSuperVector(Vec x, Vec* y);