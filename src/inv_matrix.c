#include <stdlib.h>
#include <stdio.h>

#include "inv_matrix.h"

#define INV_MATRIX_VERBOSE "MATRIX:\t\t"


Mat matrix_create_from_input(MPI_Comm comm, InputCSR* input, int verbose){

    int size;
    MPI_Comm_size(comm, &size);

    Mat Q;
    int n = input->n, m = input->m;
    MatCreate(comm, &Q);
    if(size>1){
        if(verbose) printf("%sCreating an MPI matrix on %i procs..\n", INV_MATRIX_VERBOSE, size);
        MatSetType(Q, MATMPIAIJ);
        MatSetSizes(Q, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatMPISBAIJSetPreallocation(Q, 1, 5, NULL, 2, NULL);    // must be changed if use MatMPIXXX
    } else {
        if(verbose) printf("%sCreating a SEQ matrix on %i proc..\n", INV_MATRIX_VERBOSE, size);
        MatSetType(Q, MATSEQSBAIJ);     // MATSEQAIJ?
        MatSetSizes(Q, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatSeqAIJSetPreallocation(Q, 1, input->nnz);
    }
    MatSetFromOptions(Q);
    MatSetUp(Q);
    for(int k=0; k<m; k++){
        MatSetValue(Q, input->i[k], input->j[k], input->val[k], INSERT_VALUES);
        MatSetValue(Q, input->j[k], input->i[k], input->val[k], INSERT_VALUES);
    }
    MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY);
    if(verbose) printf("%sMatrix assembled.\n", INV_MATRIX_VERBOSE);

    return Q;
}


Mat matrix_create_rw2(MPI_Comm comm, int n, double alpha, double tau, int verbose){

    int size;
    MPI_Comm_size(comm, &size);

    Mat Q;
    int Istart, Iend;
    MatCreate(comm, &Q);
    if(size>1){
        if(verbose) printf("%sCreating an MPI matrix on %i procs..\n", INV_MATRIX_VERBOSE, size);
        MatSetType(Q, MATMPIAIJ);
        MatSetSizes(Q, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatMPISBAIJSetPreallocation(Q, 1, 5, NULL, 2, NULL);
    } else {
        if(verbose) printf("%sCreating a SEQ matrix on %i proc..\n", INV_MATRIX_VERBOSE, size);
        MatSetType(Q, MATSEQSBAIJ);     // MATSEQAIJ?
        MatSetSizes(Q, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatSeqSBAIJSetPreallocation(Q, 1, 5, NULL);
    }
    MatSetFromOptions(Q);
    MatSetUp(Q);
    MatGetOwnershipRange(Q, &Istart, &Iend);
    for(int i=Istart; i<Iend; i++){
        if(i==0){
            MatSetValue(Q, i, i+0, 1 * alpha + tau, INSERT_VALUES);
            MatSetValue(Q, i, i+1, -2 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, i+2, 1 * alpha, INSERT_VALUES);
        } else if(i==1){
            MatSetValue(Q, i, i-1, -2 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, i+0, 5 * alpha + tau, INSERT_VALUES);
            MatSetValue(Q, i, i+1, -4 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, i+2, 1 * alpha, INSERT_VALUES);
        } else if(i<n-2){
            MatSetValue(Q, i, i-2, 1 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, i-1, -4 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, i+0, 6 * alpha + tau, INSERT_VALUES);
            MatSetValue(Q, i, i+1, -4 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, i+2, 1 * alpha, INSERT_VALUES);
        } else if(i<n-1){
            MatSetValue(Q, i, i-2, 1 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, i-1, -4 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, i+0, 5 * alpha + tau, INSERT_VALUES);
            MatSetValue(Q, i, i+1, -2 * alpha, INSERT_VALUES);
        } else {
            MatSetValue(Q, i, i-2, 1 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, i-1, -2 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, i+0, 1 * alpha + tau, INSERT_VALUES);
        }
    }
    MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY);
    if(verbose) printf("%sMatrix assembled.\n", INV_MATRIX_VERBOSE);

    return Q;
}


void matrix_write_binary(MPI_Comm comm, char* filename, Mat Q, int verbose){

    if(verbose) printf("%sWriting matrix..\n", INV_MATRIX_VERBOSE);
    char fullpath[256];
    snprintf(fullpath, sizeof(fullpath), "data/%s", filename);
    PetscViewer viewer;
    PetscViewerBinaryOpen(comm, fullpath, FILE_MODE_WRITE, &viewer);
    MatView(Q, viewer);
    PetscViewerDestroy(&viewer);
    if(verbose) printf("%sMatrix written.\n", INV_MATRIX_VERBOSE);
}


Mat matrix_create_submatrix(Mat Q, IS is, int verbose){

    if(verbose) printf("%sCreating a submatrix..\n", INV_MATRIX_VERBOSE);
    Mat Q_a;
    MatCreateSubMatrix(Q, is, is, MAT_INITIAL_MATRIX, &Q_a);
    if(verbose) printf("%sSubmatrix created.\n", INV_MATRIX_VERBOSE);

    return Q_a;
}