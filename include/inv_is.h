#include <petsc.h>


/**
 * @brief Load a vector from a file
 * 
 * @param comm 
 * @param filename 
 * @param is 
 * @return PetscErrorCode 
 */
PetscErrorCode ISCreateLoad(MPI_Comm comm, const char filename[], IS* is);


/**
 * @brief Create a contiguous index set from rows of a parallel matrix
 * 
 * @param comm 
 * @param A 
 * @param is 
 * @return PetscErrorCode 
 */
PetscErrorCode MatGetOwnershipISRows(MPI_Comm comm, Mat A, IS* is);


/**
 * @brief Create a block index set
 * 
 * @param comm 
 * @param is 
 * @param bs 
 * @param isb 
 * @return PetscErrorCode 
 */
PetscErrorCode ISCreateBlockIS(MPI_Comm comm, IS is, int bs, IS* isb);