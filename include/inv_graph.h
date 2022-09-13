#ifndef INV_GRAPH_H
#define INV_GRAPH_H

#include "inv_input.h"


/**
 * @brief A struct storing a graph
 * 
 */
typedef struct Graph {
    int n_vert;
    int n_edge;
    int sep_size;
    int sep_index;
    int* xadj;
    int* yadj;
    int* part;
    int* vert;
} Graph;


/**
 * @brief A constructor for Graph
 * 
 * @param n_vert 
 * @param n_edge 
 * @return Graph* 
 */
Graph* graph_create(int n_vert, int n_edge);


/**
 * @brief Create a graph from a CSR input
 * 
 * @param input 
 * @param verbose 
 * @return Graph* 
 */
Graph* graph_create_from_input(InputCSR* input, int verbose);


/**
 * @brief Create a RW2 graph of size n
 * 
 * @param n 
 * @param verbose 
 * @return Graph* 
 */
Graph* graph_create_rw2(int n, int verbose);


/**
 * @brief A destructor for Graph
 * 
 * @param graph 
 */
void graph_destroy(Graph* graph);


/**
 * @brief Print a graph
 * 
 * @param graph 
 * @param vert 
 * @param adj 
 * @param part 
 */
void graph_print(Graph* graph, int vert, int adj, int part);


/**
 * @brief Separate a graph
 * 
 * @param graph 
 * @param verbose 
 * @return Graph* 
 */
Graph* graph_separate(Graph* graph, int verbose);


/**
 * @brief Partition a graph into several parts
 * 
 * @param graph 
 * @param n_parts 
 * @param verbose 
 * @return Graph* 
 */
Graph* graph_partition(Graph* graph, int n_parts, int verbose);


/**
 * @brief Create a subgraph
 * 
 * @param graph 
 * @param idx_part 
 * @return Graph* 
 */
Graph* graph_subgraph(Graph* graph, int idx_part);


/**
 * @brief Compute the graph and its neighboring nodes
 * 
 * @param graph 
 * @param idx_part 
 * @param order 
 * @param verbose 
 * @return Graph* 
 */
Graph* graph_neighborhood(Graph* graph, int idx_part, int order, int verbose);


#endif