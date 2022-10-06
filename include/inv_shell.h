#ifndef INV_SHELL_H
#define INV_SHELL_H

#include <petsc.h>


/**
 * @brief A shell preconditioner
 * 
 */
typedef struct InvShellPC {
    Mat Q;
    PC pc;
} InvShellPC;


/**
 * @brief Create a shell preconditioner
 * 
 * @param shell 
 * @return PetscErrorCode 
 */
PetscErrorCode InvShellCreate(InvShellPC** shell);


/**
 * @brief Set up a shell preconditioner
 * 
 * @param pc 
 * @param Q 
 * @param comm 
 * @return PetscErrorCode 
 */
PetscErrorCode InvShellSetup(PC pc, Mat Q, MPI_Comm comm);


/**
 * @brief Apply a shell preconditioner to a vector
 * 
 * @param pc 
 * @param x 
 * @param y 
 * @return PetscErrorCode 
 */
PetscErrorCode InvShellApply(PC pc, Vec x, Vec y);


/**
 * @brief Apply transpose of a shell preconditioner to a vector
 * 
 * @param pc 
 * @param x 
 * @param y 
 * @return PetscErrorCode 
 */
PetscErrorCode InvShellApplyTranspose(PC pc, Vec x, Vec y);


/**
 * @brief Apply a shell preconditioner times operator to a vector
 * 
 * @param pc 
 * @param side 
 * @param x 
 * @param y 
 * @param work 
 * @return PetscErrorCode 
 */
PetscErrorCode InvShellApplyBA(PC pc, PCSide side, Vec x, Vec y, Vec work);


/**
 * @brief Apply a shell preconditioner symmetrically from the left (does nothing for sampler KSP)
 * 
 * @param pc 
 * @param x 
 * @param y 
 * @return PetscErrorCode 
 */
PetscErrorCode InvShellApplyLeft(PC pc, Vec x, Vec y);


/**
 * @brief Apply a shell preconditioner symmetrically from the right
 * 
 * @param pc 
 * @param x 
 * @param y 
 * @return PetscErrorCode 
 */
PetscErrorCode InvShellApplyRight(PC pc, Vec x, Vec y);


/**
 * @brief Destroy a shell preconditioner
 * 
 * @param pc 
 * @return PetscErrorCode 
 */
PetscErrorCode InvShellDestroy(PC pc);


#endif