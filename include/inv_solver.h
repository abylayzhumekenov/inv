#ifndef INV_SOLVER_H
#define INV_SOLVER_H

#include <petsc.h>

#define INV_SOLVER_VERBOSE "SOLVER:\t\t"


/**
 * @brief Create a KSP solver
 * 
 * @param comm 
 * @param Q 
 * @param max_niter 
 * @param verbose 
 * @return KSP 
 */
KSP solver_ksp(MPI_Comm comm, Mat Q, int max_niter, int verbose);


/**
 * @brief Compute the sample contribution
 * 
 * @param ksp 
 * @param Q 
 * @param x 
 * @param is_x 
 * @param is_a 
 * @param n_samples 
 * @param verbose 
 * @return Vec 
 */
Vec solver_sample_contribution(KSP ksp, Mat Q, Vec x, IS is_x, IS is_a, int verbose);


/**
 * @brief Solve the base case
 * 
 * @param comm 
 * @param Q 
 * @param verbose 
 * @return Vec 
 */
Vec solver_base_case(MPI_Comm comm, Mat Q, int verbose);


/**
 * @brief Assemble the solution
 * 
 * @param d 
 * @param is_part 
 * @param is_a 
 * @param n 
 * @param verbose 
 * @return Vec 
 */
Vec solver_assemble_solution(Vec d, IS is_part, IS is_a, int n, int verbose);

#endif