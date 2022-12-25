#ifndef INV_GRAPH_H
#define INV_GRAPH_H


typedef struct InvGraph {
    int n_vert;
    int n_edge;
    int n_part;
    int i_part;
    int* xadj;
    int* yadj;
    int* part;
    int* vert;
    int* offset;
} InvGraph;


InvGraph* InvGraphCreate(char* filename, int verbose);


void InvGraphDestroy(InvGraph* graph, int verbose);


void InvGraphPrint(InvGraph* graph, int verbose);


void InvGraphPartition(InvGraph* graph, int n_part, int verbose);


void InvGraphNeighbor(InvGraph* graph, int n_neighbor, int verbose);


#endif