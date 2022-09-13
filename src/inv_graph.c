#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <petsc.h>
#include <metis.h>

#include "inv_graph.h"

#define INV_GRAPH_VERBOSE "GRAPH:\t\t"


Graph* graph_create(int n_vert, int n_edge){

    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->n_vert = n_vert;
    graph->n_edge = n_edge;
    graph->sep_index = 0;
    graph->sep_size = 0;
    graph->xadj = (int*)malloc((n_vert+1) * sizeof(int));
    graph->yadj = (int*)malloc(n_edge * sizeof(int));
    graph->vert = (int*)malloc(n_vert * sizeof(int));
    graph->part = (int*)calloc(n_vert, sizeof(int));

    for(int i=0; i<n_vert; i++) graph->vert[i] = i;

    return graph;
}


Graph* graph_create_from_input(InputCSR* input, int verbose){

    if(verbose) printf("%sCreating a graph..\n", INV_GRAPH_VERBOSE);
    int n_vert = input->n;
    int n_edge = 2 * (input->m - input->n);
    Graph* graph = graph_create(n_vert, n_edge);
    memcpy(graph->xadj, input->xadj, (n_vert+1) * sizeof(int));
    memcpy(graph->yadj, input->yadj, n_edge * sizeof(int));

    return graph;
}


Graph* graph_create_rw2(int n, int verbose){

    if(verbose) printf("%sCreating a graph..\n", INV_GRAPH_VERBOSE);
    int n_vert = n, n_edge = 4*n-6;
    Graph* graph = graph_create(n_vert, n_edge);
    for(int i=0; i<n_vert+1; i++){
        if(i==0){
            graph->xadj[i] = 0;
        } else if(i==1){
            graph->xadj[i] = graph->xadj[i-1] + 2;
            graph->yadj[graph->xadj[i-1]+0] = i-1+1;
            graph->yadj[graph->xadj[i-1]+1] = i-1+2;
        } else if(i==2){
            graph->xadj[i] = graph->xadj[i-1] + 3;
            graph->yadj[graph->xadj[i-1]+0] = i-1-1;
            graph->yadj[graph->xadj[i-1]+1] = i-1+1;
            graph->yadj[graph->xadj[i-1]+2] = i-1+2;
        } else if(i<n_vert-1){
            graph->xadj[i] = graph->xadj[i-1] + 4;
            graph->yadj[graph->xadj[i-1]+0] = i-1-2;
            graph->yadj[graph->xadj[i-1]+1] = i-1-1;
            graph->yadj[graph->xadj[i-1]+2] = i-1+1;
            graph->yadj[graph->xadj[i-1]+3] = i-1+2;
        } else if(i==n_vert-1){
            graph->xadj[i] = graph->xadj[i-1] + 3;
            graph->yadj[graph->xadj[i-1]+0] = i-1-2;
            graph->yadj[graph->xadj[i-1]+1] = i-1-1;
            graph->yadj[graph->xadj[i-1]+2] = i-1+1;
        } else {
            graph->xadj[i] = graph->xadj[i-1] + 2;
            graph->yadj[graph->xadj[i-1]+0] = i-1-2;
            graph->yadj[graph->xadj[i-1]+1] = i-1-1;
        }
    }

    return graph;
}


void graph_destroy(Graph* graph){

    free(graph->xadj);
    free(graph->yadj);
    free(graph->vert);
    free(graph->part);
    free(graph);

    return;
}


void graph_print(Graph* graph, int vert, int adj, int part){

    if(vert){
        for(int i=0; i<graph->n_vert; i++){
            if(i%20==0) printf("\n");
            printf("%i\t", graph->vert[i]);
        }
        printf("\n");
    }
    if(adj){
        for(int i=0; i<graph->n_vert+1; i++){
            if(i%20==0) printf("\n");
            printf("%i\t", graph->xadj[i]);
        }
        printf("\n");
        for(int i=0; i<graph->n_edge; i++){
            if(i%20==0) printf("\n");
            printf("%i\t", graph->yadj[i]);
        }
        printf("\n");
    }
    if(part){
        for(int i=0; i<graph->n_vert; i++){
            if(i%20==0) printf("\n");
            printf("%i\t", graph->part[i]);
        }
        printf("\n");
    }
}


Graph* graph_separate(Graph* graph, int verbose){

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size > 1){
        if(verbose) printf("%sSeparating a graph..\n", INV_GRAPH_VERBOSE);
        METIS_ComputeVertexSeparator(&graph->n_vert, graph->xadj, graph->yadj, NULL, NULL, &graph->sep_size, graph->part);
        if(verbose) printf("%sSeparation finished.\n", INV_GRAPH_VERBOSE);
    } else {
        if(verbose) printf("%sNothing to separate.\n", INV_GRAPH_VERBOSE);
    }

    return graph;
}


Graph* graph_partition(Graph* graph, idx_t n_parts, int verbose){

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size > 1){
        idx_t n_con = 1, objval = 0;
        if(verbose) printf("%sPartitioning a graph..\n", INV_GRAPH_VERBOSE);
        METIS_PartGraphRecursive(&graph->n_vert, &n_con, graph->xadj, graph->yadj, NULL, NULL, NULL, &n_parts, NULL, NULL, NULL, &objval, graph->part);
        if(verbose) printf("%sPartitioning finished.\n", INV_GRAPH_VERBOSE);
    } else {
        if(verbose) printf("%sNothing to partition.\n", INV_GRAPH_VERBOSE);
    }

    return graph;
}


Graph* graph_neighborhood(Graph* graph, int idx_part, int order, int verbose){

    if(verbose) printf("%sComputing neighbors..\n", INV_GRAPH_VERBOSE);
    for(int p=0; p<order; p++){
        int part_temp[graph->n_vert];
        for(int i=0; i<graph->n_vert; i++) part_temp[i] = graph->part[i];
        for(int i=0; i<graph->n_vert; i++){
            if(graph->part[i]==idx_part){
                for(int j=graph->xadj[i]; j<graph->xadj[i+1]; j++){
                    part_temp[graph->yadj[j]] = idx_part;
                }
            }
        }
        for(int i=0; i<graph->n_vert; i++) graph->part[i] = part_temp[i];
    }
    if(verbose) printf("%sNeighborhood computed.\n", INV_GRAPH_VERBOSE);
    
    return graph;
}
