#ifndef INV_IS_H
#define INV_IS_H

#include <petsc.h>
#include "inv_graph.h"


/**
 * @brief Index set mapping
 * 
 */
typedef struct InvIS {
    int n_vert;
    int n_part;
    int* offset;
    int* n_local;
    int* index;
    int* part;
} InvIS;


/**
 * @brief Create a mapping
 * 
 * @param n_vert 
 * @param n_part 
 * @return InvIS* 
 */
InvIS* InvISCreate(int n_vert, int n_part);


/**
 * @brief Create a mapping from an initial graph to a reordered one
 * 
 * @param graph 
 * @param verbose 
 * @return InvIS* 
 */
InvIS* InvISCreateInitialToGlobal(InvGraph* graph, int verbose);


/**
 * @brief Create a mapping from a reordered graph to an initial one
 * 
 * @param graph 
 * @param verbose 
 * @return InvIS* 
 */
InvIS* InvISCreateGlobalToInitial(InvGraph* graph, int verbose);


/**
 * @brief Invert a mapping (works only for Initial-Global reordering)
 * 
 * @param mapping 
 * @param verbose 
 * @return InvIS* 
 */
InvIS* InvISCreateInverse(InvIS* mapping, int verbose);


/**
 * @brief Create a mapping from an initial graph to a reordered local subgraph
 * 
 * @param graph 
 * @param verbose 
 * @return InvIS* 
 */
InvIS* InvISCreateInitialToLocal(InvGraph* graph, int verbose);


/**
 * @brief Create a mapping from a reordered local subgraph to an initial graph (breaks for nonlocal indices)
 * 
 * @param graph 
 * @param i_part 
 * @param verbose 
 * @return InvIS* 
 */
InvIS* InvISCreateLocalToInitial(InvGraph* graph, int i_part, int verbose);


/**
 * @brief Destroy a mapping
 * 
 * @param mapping 
 * @param verbose 
 */
void InvISDestroy(InvIS* mapping, int verbose);


/**
 * @brief print a mapping
 * 
 * @param mapping 
 * @param verbose 
 */
void InvISPrint(InvIS* mapping, int verbose);


#endif