#ifndef INV_IS_H
#define INV_IS_H

#include <petsc.h>
#include "inv_graph.h"


/**
 * @brief Create an index set from a partition
 * 
 * @param comm 
 * @param graph 
 * @param idx_part 
 * @param verbose 
 * @return IS 
 */
IS is_create_from_partition(MPI_Comm comm, Graph* graph, int idx_part, int verbose);

#endif