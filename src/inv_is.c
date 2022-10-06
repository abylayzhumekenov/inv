#include "inv_is.h"

#define INV_IS_VERBOSE "IS:\t\t"


InvIS* InvISCreate(int n_vert, int n_part){

    InvIS* mapping = (InvIS*)malloc(sizeof(InvIS));
    mapping->n_vert = n_vert;
    mapping->n_part = n_part;
    mapping->offset = (int*)calloc(n_part, sizeof(int));
    mapping->n_local = (int*)calloc(n_part, sizeof(int));
    mapping->index = (int*)calloc(n_vert, sizeof(int));
    mapping->part = (int*)calloc(n_vert, sizeof(int));

    return mapping;
}


InvIS* InvISCreateInitialToGlobal(InvGraph* graph, int verbose){

    int n_vert = graph->n_vert, n_part = graph->n_part;
    InvIS* mapping = InvISCreate(n_vert, n_part);

    int temp[n_part];
    for(int i=0; i<n_part; i++) temp[i] = 0;
    for(int i=0; i<n_vert; i++) mapping->n_local[graph->part[i]]++;
    for(int i=1; i<n_part; i++) mapping->offset[i] = mapping->offset[i-1] + mapping->n_local[i-1];
    for(int i=0; i<n_vert; i++){
        mapping->part[i] = graph->part[i];
        mapping->index[i] = mapping->offset[graph->part[i]] + temp[graph->part[i]];
        temp[graph->part[i]]++;
    }
    if(verbose) printf("%sMapping created.\n", INV_IS_VERBOSE);

    return mapping;
}

InvIS* InvISCreateGlobalToInitial(InvGraph* graph, int verbose){

    int n_vert = graph->n_vert, n_part = graph->n_part;
    InvIS* mapping = InvISCreate(n_vert, n_part);

    int temp[n_part];
    for(int i=0; i<n_part; i++) temp[i] = 0;
    for(int i=0; i<n_vert; i++) mapping->n_local[graph->part[i]]++;
    for(int i=1; i<n_part; i++) mapping->offset[i] = mapping->offset[i-1] + mapping->n_local[i-1];
    for(int i=0; i<n_vert; i++){
        mapping->part[i] = graph->part[i];
        mapping->index[mapping->offset[graph->part[i]] + temp[graph->part[i]]] = i;
        temp[graph->part[i]]++;
    }
    if(verbose) printf("%sMapping created.\n", INV_IS_VERBOSE);

    return mapping;
}


InvIS* InvISCreateInverse(InvIS* mapping, int verbose){

    InvIS* inverse = InvISCreate(mapping->n_vert, mapping->n_part);
    for(int i=0; i<mapping->n_vert; i++){
        inverse->index[mapping->index[i]] = i;
    }
    if(verbose) printf("%sMapping inverted.\n", INV_IS_VERBOSE);

    return inverse;
}


InvIS* InvISCreateInitialToLocal(InvGraph* graph, int verbose){

    int n_vert = graph->n_vert, n_part = graph->n_part;
    InvIS* mapping = InvISCreate(n_vert, n_part);

    int temp[n_part];
    for(int i=0; i<n_part; i++) temp[i] = 0;
    for(int i=0; i<n_vert; i++) mapping->n_local[graph->part[i]]++;
    for(int i=1; i<n_part; i++) mapping->offset[i] = mapping->offset[i-1] + mapping->n_local[i-1];
    for(int i=0; i<n_vert; i++){
        mapping->part[i] = graph->part[i];
        mapping->index[i] = temp[graph->part[i]];
        temp[graph->part[i]]++;
    }
    if(verbose) printf("%sMapping created.\n", INV_IS_VERBOSE);

    return mapping;
}


InvIS* InvISCreateLocalToInitial(InvGraph* graph, int i_part, int verbose){

    int n_vert = graph->n_vert, n_part = graph->n_part;
    InvIS* mapping = InvISCreate(n_vert, n_part);

    int temp[n_part];
    for(int i=0; i<n_part; i++) temp[i] = 0;
    for(int i=0; i<n_vert; i++) mapping->n_local[graph->part[i]]++;
    for(int i=1; i<n_part; i++) mapping->offset[i] = mapping->offset[i-1] + mapping->n_local[i-1];
    for(int i=0; i<n_vert; i++){
        mapping->part[i] = graph->part[i];
        if(graph->part[i]==i_part){
            mapping->index[temp[graph->part[i]]] = i;
            temp[graph->part[i]]++;
        }
    }
    if(verbose) printf("%sMapping created.\n", INV_IS_VERBOSE);

    return mapping;
}


void InvISDestroy(InvIS* mapping, int verbose){

    free(mapping->index);
    free(mapping);
    if(verbose) printf("%sMapping destroyed.\n", INV_IS_VERBOSE);
}


void InvISPrint(InvIS* mapping, int verbose){

    if(verbose){
        printf("Mapping:\n");
        for(int i=0; i<mapping->n_vert; i++){
            if(i<10){
                printf("%i\t", mapping->index[i]);
            } else if(i==10){
                printf("...");
                break;
            }
        }
        printf("\n");
    }
}