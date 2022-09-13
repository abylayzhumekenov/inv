#include "inv_is.h"

#define INV_IS_VERBOSE "IS:\t\t"


IS is_create_from_partition(MPI_Comm comm, Graph* graph, int idx_part, int verbose){

    if(verbose) printf("%sCreating an index set..\n", INV_IS_VERBOSE);
    int n_is = 0;
    for(int i=0; i<graph->n_vert; i++){
        if(graph->part[i]==idx_part) n_is++;
    }
    int is_temp[n_is], k = 0;
    for(int i=0; i<graph->n_vert; i++){
        if(graph->part[i]==idx_part) {
            is_temp[k] = i; 
            k++;
        }
    }
    IS is;
    ISCreateGeneral(comm, n_is, is_temp, PETSC_COPY_VALUES, &is);
    if(verbose) printf("%sIndex set created.\n", INV_IS_VERBOSE);

    return is;
}