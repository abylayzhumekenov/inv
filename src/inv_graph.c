#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>
#include <metis.h>

#include "inv_graph.h"

#define INV_GRAPH_VERBOSE "GRAPH:\t\t"


InvGraph* InvGraphCreate(int n_vert, int n_edge, int verbose){

    InvGraph *graph = (InvGraph*)malloc(sizeof(InvGraph));
    graph->n_vert = n_vert;
    graph->n_edge = n_edge;
    graph->n_part = 1;
    graph->xadj = (int*)calloc(n_vert+1, sizeof(int));
    graph->yadj = (int*)calloc(n_edge, sizeof(int));
    graph->vert = (int*)calloc(n_vert, sizeof(int));
    graph->part = (int*)calloc(n_vert, sizeof(int));
    for(int i=0; i<n_vert; i++) graph->vert[i] = i;

    return graph;
}


InvGraph* InvGraphCreateFromFile(char* filename, int verbose){

    /* Open a file */
    char fullpath[256];
    snprintf(fullpath, sizeof(fullpath), "data/%s", filename);
    FILE* file = fopen(fullpath, "r");
    if(!file){
        if(verbose) printf("%sFile not found.\n", INV_GRAPH_VERBOSE);
        fclose(file);
        return NULL;
    } else {
        if(verbose) printf("%sReading a file \"%s\"...\n", INV_GRAPH_VERBOSE, fullpath);
    }

    /* Read the header */
    char* line = NULL;
    size_t line_len;
    int n, m, i, j;
    double val;
    getline(&line, &line_len, file);
    fscanf(file, "%i%i%i\n", &n, &n, &m);
    if(verbose) printf("%s%s", INV_GRAPH_VERBOSE, line+2);
    if(verbose) printf("%s%i rows, %i non-zeros.\n", INV_GRAPH_VERBOSE, n, m);

    /* Compute number of nonzeros */
    int nnz[n], nnz_temp[n];
    for(int k=0; k<n; k++){
        nnz[k] = nnz_temp[k] = 0;
    }
    for(int k=0; k<m; k++){
        fscanf(file, "%i%i%lf\n", &i, &j, &val);
        nnz[j-1] += (i!=j);
        nnz[i-1] += (i!=j);
    }

    /* Create a graph */
    if(verbose) printf("%sCreating a graph...\n", INV_GRAPH_VERBOSE);
    int n_vert = n, n_edge = 2*(m-n);
    InvGraph *graph = InvGraphCreate(n_vert, n_edge, verbose);

    /* Fill adjacency lists */
    fseek(file, 0, SEEK_SET);
    getline(&line, &line_len, file);
    fscanf(file, "%i%i%i\n", &n, &n, &m);
    for(int k=0; k<n; k++){
        graph->xadj[k+1] = graph->xadj[k] + nnz[k];
    }
    for(int k=0; k<m; k++){
        fscanf(file, "%i%i%lf\n", &i, &j, &val);
        graph->yadj[graph->xadj[i-1] + nnz_temp[i-1]] += (i!=j) * (j-1);
        graph->yadj[graph->xadj[j-1] + nnz_temp[j-1]] += (i!=j) * (i-1);
        nnz_temp[j-1] += (i!=j);
        nnz_temp[i-1] += (i!=j);
    }
    fclose(file);
    if(verbose) printf("%sGraph created.\n", INV_GRAPH_VERBOSE);

    return graph;
}


InvGraph* InvGraphCreateRW2(int n_vert, int cyclic, int verbose){

    /* Create a graph */
    if(verbose) printf("%sCreating a graph...\n", INV_GRAPH_VERBOSE);
    int n_edge = n_vert*4 - (cyclic==0)*6;
    InvGraph *graph = InvGraphCreate(n_vert, n_edge, verbose);

    /* Fill adjacency lists */
    for(int k=0; k<n_vert; k++){
        if(k==0){
            graph->xadj[k+1] = graph->xadj[k] + 2;
            graph->yadj[graph->xadj[k]+0] = k+1;
            graph->yadj[graph->xadj[k]+1] = k+2;
            if(cyclic){
                graph->xadj[k+1] += 2;
                graph->yadj[graph->xadj[k]+2] = (k-2+n_vert) % n_vert;
                graph->yadj[graph->xadj[k]+3] = (k-1+n_vert) % n_vert;
            }
        } else if(k==1){
            graph->xadj[k+1] = graph->xadj[k] + 3;
            graph->yadj[graph->xadj[k]+0] = k-1;
            graph->yadj[graph->xadj[k]+1] = k+1;
            graph->yadj[graph->xadj[k]+2] = k+2;
            if(cyclic){
                graph->xadj[k+1] += 1;
                graph->yadj[graph->xadj[k]+3] = (k-2+n_vert) % n_vert;
            }
        } else if(k<n_vert-2){
            graph->xadj[k+1] = graph->xadj[k] + 4;
            graph->yadj[graph->xadj[k]+0] = k-2;
            graph->yadj[graph->xadj[k]+1] = k-1;
            graph->yadj[graph->xadj[k]+2] = k+1;
            graph->yadj[graph->xadj[k]+3] = k+2;
        } else if(k<n_vert-1){
            graph->xadj[k+1] = graph->xadj[k] + 3;
            graph->yadj[graph->xadj[k]+0] = k-2;
            graph->yadj[graph->xadj[k]+1] = k-1;
            graph->yadj[graph->xadj[k]+2] = k+1;
            if(cyclic){
                graph->xadj[k+1] += 1;
                graph->yadj[graph->xadj[k]+3] = (k+2+n_vert) % n_vert;
            }
        } else if(k<n_vert){
            graph->xadj[k+1] = graph->xadj[k] + 2;
            graph->yadj[graph->xadj[k]+0] = k-2;
            graph->yadj[graph->xadj[k]+1] = k-1;
            if(cyclic){
                graph->xadj[k+1] += 2;
                graph->yadj[graph->xadj[k]+2] = (k+1+n_vert) % n_vert;
                graph->yadj[graph->xadj[k]+3] = (k+2+n_vert) % n_vert;
            }
        }
    }
    if(verbose) printf("%sGraph created.\n", INV_GRAPH_VERBOSE);

    return graph;
}


void InvGraphDestroy(InvGraph* graph, int verbose){

    free(graph->xadj);
    free(graph->yadj);
    free(graph->vert);
    free(graph->part);
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


InvGraph* InvGraphPartition(InvGraph* graph, int n_part, int verbose){

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(n_part == 1){
        if(verbose) printf("%sPartitioning aborted. Only one partition.\n", INV_GRAPH_VERBOSE);
    } else if(n_part <= size){
        int n_con = 1, objval = 0;
        if(verbose) printf("%sPartitioning a graph..\n", INV_GRAPH_VERBOSE);
        METIS_PartGraphRecursive(&graph->n_vert, &n_con, graph->xadj, graph->yadj, NULL, NULL, NULL, &n_part, NULL, NULL, NULL, &objval, graph->part);
        // METIS_PartGraphKway(&graph->n_vert, &n_con, graph->xadj, graph->yadj, NULL, NULL, NULL, &n_part, NULL, NULL, NULL, &objval, graph->part);
        graph->n_part = n_part;
        if(verbose) printf("%sPartitioning finished.\n", INV_GRAPH_VERBOSE);
    } else {
        if(verbose) printf("%sPartitioning aborted. World size is less than number of partitions.\n", INV_GRAPH_VERBOSE);
    }

    return graph;
}


InvGraph* InvGraphNeighborhood(InvGraph* graph, int i_part, int n_neighbor, int verbose){

    if(verbose) printf("%sComputing neighbors...\n", INV_GRAPH_VERBOSE);
    for(int p=0; p<n_neighbor; p++){
        int part_temp[graph->n_vert];
        for(int i=0; i<graph->n_vert; i++) part_temp[i] = graph->part[i];
        for(int i=0; i<graph->n_vert; i++){
            if(graph->part[i]==i_part){
                for(int j=graph->xadj[i]; j<graph->xadj[i+1]; j++){
                    part_temp[graph->yadj[j]] = i_part;
                }
            }
        }
        for(int i=0; i<graph->n_vert; i++) graph->part[i] = part_temp[i];
    }
    if(verbose) printf("%sNeighborhood computed.\n", INV_GRAPH_VERBOSE);

    return graph;
}
