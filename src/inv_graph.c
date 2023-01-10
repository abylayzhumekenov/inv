#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>
#include <metis.h>

#include "inv_graph.h"


#define INV_GRAPH_VERBOSE "GRAPH:\t\t"

uint32_t htobe32(uint32_t x) {
	union { uint32_t u32; uint8_t v[4]; } ret;
	ret.v[0] = (uint8_t)(x >> 24);
	ret.v[1] = (uint8_t)(x >> 16);
	ret.v[2] = (uint8_t)(x >> 8);
	ret.v[3] = (uint8_t) x;
	return ret.u32;
}


InvGraph* InvGraphCreate(char* filename, int verbose){

    /* Open a file */
    if(verbose) printf("%sOpening a file \"%s\"...\n", INV_GRAPH_VERBOSE, filename);
    FILE* file = fopen(filename, "rb");
    if(!file){
        if(verbose) printf("%sFile not found.\n", INV_GRAPH_VERBOSE);
        fclose(file);
        return NULL;
    }

    /* Read the header */
    if(verbose) printf("%sReading the header...\n", INV_GRAPH_VERBOSE);
    unsigned int graph_header[4];
    fread(graph_header, 4, 4, file);
    int n = htobe32(graph_header[2]), nnz = htobe32(graph_header[3]);
    int n_vert = n, n_edge = 2*(nnz-n);

    /* Allocate memory */
    if(verbose) printf("%sAllocating memory...\n", INV_GRAPH_VERBOSE);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    InvGraph* graph;
    graph = (InvGraph*)malloc(sizeof(InvGraph));
    graph->n_vert = n_vert;
    graph->n_edge = n_edge;
    graph->n_part = 1;
    graph->i_part = rank;
    graph->xadj = (int*)calloc(n_vert+1, sizeof(int));
    graph->yadj = (int*)calloc(n_edge, sizeof(int));
    graph->part = (int*)calloc(n_vert, sizeof(int));
    graph->vert = (int*)calloc(n_vert, sizeof(int));
    graph->offset = (int*)calloc(4, sizeof(int));

    /* Fill the adjacency */
    if(verbose) printf("%sFilling the graph...\n", INV_GRAPH_VERBOSE);
    unsigned int buf_nnz_row[n], buf_nnz_i[nnz];
    fread(buf_nnz_row, n, 4, file);
    fread(buf_nnz_i, nnz, 4, file);
    int k = 0, l = 0;
    for(int i=0; i<n; i++){
        graph->xadj[i+1] = graph->xadj[i] + htobe32(buf_nnz_row[i]) - 1;
        for(int j=0; j<htobe32(buf_nnz_row[i]); j++){
            if(i != htobe32(buf_nnz_i[k])){
                graph->yadj[l] = htobe32(buf_nnz_i[k]);
                l++;
            }
            k++;
        }
    }

    fclose(file);
    if(verbose) printf("%sGraph created.\n", INV_GRAPH_VERBOSE);

    return graph;
}


void InvGraphDestroy(InvGraph* graph, int verbose){

    free(graph->xadj);
    free(graph->yadj);
    free(graph->part);
    free(graph->vert);
    free(graph->offset);
    free(graph);
    if(verbose) printf("%sGraph destroyed.\n", INV_GRAPH_VERBOSE);
}


void InvGraphPrint(InvGraph* graph, int verbose){

    if(verbose){
        printf("Graph:\n");
        for(int i=0; i<graph->n_vert; i++){
            if(i<5 || i>graph->n_vert-6){
                printf("%i\t[%i\t,%i\t):\t", i, graph->xadj[i], graph->xadj[i+1]);
                for(int j=graph->xadj[i]; j<graph->xadj[i+1]; j++){
                    printf("%i\t", graph->yadj[j]);
                }
                printf("\n");
            } else if(i==5){
                printf("...\n...\n...\n");
            }
        }
        printf("Partition:\n");
        for(int k=0; k<graph->n_part; k++){
            printf("%i\t:\t", k);
            int j=0;
            for(int i=0; i<graph->n_vert; i++){
                if(graph->part[i]==k && j<10){
                    printf("%i\t", i);
                    j++;
                } else if(graph->part[i]==k && j==10){
                    printf("...");
                    j++;
                    break;
                }
            }
            printf("\n");
        }
    }
}


void InvGraphPartition(InvGraph* graph, int n_part, int verbose){

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(n_part == 1){
        if(verbose) printf("%sPartitioning aborted. Only one partition.\n", INV_GRAPH_VERBOSE);
    } else if(n_part <= size){
        int n_con = 1, objval = 0;
        if(verbose) printf("%sPartitioning a graph..\n", INV_GRAPH_VERBOSE);
        // METIS_PartGraphRecursive(&graph->n_vert, &n_con, graph->xadj, graph->yadj, NULL, NULL, NULL, &n_part, NULL, NULL, NULL, &objval, graph->part);
        METIS_PartGraphKway(&graph->n_vert, &n_con, graph->xadj, graph->yadj, NULL, NULL, NULL, &n_part, NULL, NULL, NULL, &objval, graph->part);
        int j = 0;
        for(int i=0; i<graph->n_vert; i++){
            if(graph->part[i]==rank){
                graph->vert[j] = i;
                j++;
            }
        }
        graph->offset[1] = j;
        graph->offset[2] = 0;
        graph->offset[3] = j;
        graph->n_part = n_part;
        if(verbose) printf("%sPartitioning finished.\n", INV_GRAPH_VERBOSE);
    } else {
        if(verbose) printf("%sPartitioning aborted. World size is less than number of partitions.\n", INV_GRAPH_VERBOSE);
    }
}


void InvGraphNeighbor(InvGraph* graph, int n_neighbor, int verbose){

    if(verbose) printf("%sComputing neighbors...\n", INV_GRAPH_VERBOSE);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for(int p=0; p<n_neighbor+1; p++){
        int k = 0;
        int part_temp[graph->n_vert];
        for(int i=0; i<graph->n_vert; i++) part_temp[i] = graph->part[i];
        for(int i=graph->offset[2]; i<graph->offset[3]; i++){
            for(int j=graph->xadj[graph->vert[i]]; j<graph->xadj[graph->vert[i]+1]; j++){
                if(graph->part[graph->yadj[j]]!=rank){
                    graph->vert[graph->offset[3]+k] = graph->yadj[j];
                    k += part_temp[graph->yadj[j]]!=rank;
                    part_temp[graph->yadj[j]] = rank;
                }
            }
        }
        for(int i=0; i<graph->n_vert; i++) graph->part[i] = part_temp[i];
        graph->offset[2] = graph->offset[3];
        graph->offset[3] += k;
    }
    if(verbose) printf("%sNeighborhood computed.\n", INV_GRAPH_VERBOSE);
}
