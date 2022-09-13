#include <stdlib.h>
#include <stdio.h>

#include "inv_input.h"

#define INV_INPUT_VERBOSE "INPUT:\t\t"


InputCSR* input_create(int n, int m){

    InputCSR* input = (InputCSR*)malloc(sizeof(InputCSR));
    input->n = n;
    input->m = m;
    input->i = (int*)malloc(m * sizeof(int));
    input->j = (int*)malloc(m * sizeof(int));
    input->val = (double*)malloc(m * sizeof(double));
    input->nnz = (int*)calloc(n, sizeof(int));
    input->xadj = (int*)calloc(n+1, sizeof(int));
    input->yadj = (int*)malloc(2*(m-n) * sizeof(int));

    return input;
}


void input_destroy(InputCSR* input){

    free(input->i);
    free(input->j);
    free(input->val);
    free(input->nnz);
    free(input->xadj);
    free(input->yadj);
    free(input);

    return;
}


InputCSR* input_read(char* filename, int verbose){

    char fullpath[256];
    snprintf(fullpath, sizeof(fullpath), "data/%s", filename);
    FILE* file;
    file = fopen(fullpath, "r");
    if(!file) {
        printf("%sFile not found.\n", INV_INPUT_VERBOSE);
        fclose(file);
        return NULL;
    } else if(verbose) {
        printf("%sReading a file \"%s\"..\n", INV_INPUT_VERBOSE, fullpath);
    }

    char* line = NULL;
    size_t line_len;
    int n, m, i, j;
    double val;
    getline(&line, &line_len, file);
    fscanf(file, "%i%i%i\n", &n, &n, &m);
    if(verbose) printf("%s%s", INV_INPUT_VERBOSE, line+2);
    if(verbose) printf("%s%i rows, %i non-zeros..\n", INV_INPUT_VERBOSE, n, m);
    
    if(verbose) printf("%sAllocating memory..\n", INV_INPUT_VERBOSE);
    InputCSR* input = input_create(n, m);
    for(int k=0; k<m; k++){
        fscanf(file, "%i%i%lf\n", &i, &j, &val);
        input->i[k] = i - 1;
        input->j[k] = j - 1;
        input->val[k] = val;
        input->nnz[i-1]++;
        if(i!=j) input->nnz[j-1]++;
    }
    if(verbose) printf("%sClosing a file..\n", INV_INPUT_VERBOSE);
    fclose(file);

    if(verbose) printf("%sCreating adjacency lists..\n", INV_INPUT_VERBOSE);
    int nnz_temp[n];
    for(int k=0; k<n; k++){
        input->xadj[k+1] = input->xadj[k] + input->nnz[k] - 1;
        nnz_temp[k] = 0;
    }
    for(int k=0; k<m; k++){
        if(input->i[k]!=input->j[k]){
            input->yadj[input->xadj[input->i[k]] + nnz_temp[input->i[k]]] = input->j[k];
            input->yadj[input->xadj[input->j[k]] + nnz_temp[input->j[k]]] = input->i[k];
            nnz_temp[input->i[k]]++;
            nnz_temp[input->j[k]]++;
        }
    }
    if(verbose) printf("%sSuccessfully finished.\n", INV_INPUT_VERBOSE);

    return input;
}
