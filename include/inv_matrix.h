#ifndef INV_MATRIX_H
#define INV_MATRIX_H

#include <petsc.h>
#include "inv_input.h"
#include "inv_graph.h"
#include "inv_is.h"


/**
 * @brief Create a matrix from a CSR input
 * 
 * @param comm 
 * @param input 
 * @param verbose 
 * @return Mat 
 */
Mat matrix_create_from_input(MPI_Comm comm, InputCSR* input, int verbose);


/**
 * @brief Create a RW2 matrix of size n
 * 
 * @param comm 
 * @param n 
 * @param alpha 
 * @param tau 
 * @param verbose 
 * @return Mat 
 */
Mat matrix_create_rw2(MPI_Comm comm, int n, double alpha, double tau, int verbose);


/**
 * @brief Write a matrix to a binary file
 * 
 * @param comm 
 * @param filename 
 * @param Q 
 * @param verbose 
 */
void matrix_write_binary(MPI_Comm comm, char* filename, Mat Q, int verbose);


/**
 * @brief Create a submatrix from a partitioning
 * 
 * @param comm 
 * @param Q 
 * @param graph 
 * @param idx_part 
 * @param verbose 
 * @return Mat 
 */
Mat matrix_create_submatrix(Mat Q, IS is, int verbose);


#endif