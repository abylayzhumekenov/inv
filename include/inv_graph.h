#ifndef INV_GRAPH_H
#define INV_GRAPH_H


/**
 * @brief A graph structure
 * 
 */
typedef struct InvGraph {
    int n_vert;
    int n_edge;
    int n_part;
    int* xadj;
    int* yadj;
    int* part;
    int* vert;
} InvGraph;


/**
 * @brief Create a graph
 * 
 * @param n_vert 
 * @param n_edge 
 * @param verbose 
 * @return InvGraph* 
 */
InvGraph* InvGraphCreate(int n_vert, int n_edge, int verbose);


/**
 * @brief Create a graph from an input file in a MatrixMarket format
 * 
 * @param filename 
 * @param verbose 
 * @return InvGraph* 
 */
InvGraph* InvGraphCreateFromFile(char* filename, int verbose);


/**
 * @brief Create a graph for a RW2 model
 * 
 * @param n_vert 
 * @param cyclic 
 * @param verbose 
 * @return InvGraph* 
 */
InvGraph* InvGraphCreateRW2(int n_vert, int cyclic, int verbose);


/**
 * @brief Destroy a graph
 * 
 * @param graph 
 * @param verbose 
 */
void InvGraphDestroy(InvGraph* graph, int verbose);


/**
 * @brief Print a graph
 * 
 * @param graph 
 * @param verbose 
 */
void InvGraphPrint(InvGraph* graph, int verbose);


/**
 * @brief Partition a graph into several parts
 * 
 * @param graph 
 * @param n_part 
 * @param verbose 
 * @return InvGraph* 
 */
InvGraph* InvGraphPartition(InvGraph* graph, int n_part, int verbose);


/**
 * @brief Compute a neighborhood of a graph partition
 * 
 * @param graph 
 * @param i_part 
 * @param n_neighbor 
 * @param verbose 
 * @return InvGraph* 
 */
InvGraph* InvGraphNeighborhood(InvGraph* graph, int i_part, int n_neighbor, int verbose);


#endif