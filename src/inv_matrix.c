#include "inv_matrix.h"

#define INV_MATRIX_VERBOSE "MATRIX:\t\t"


Mat InvMatrixCreateFromFile(MPI_Comm comm, char* filename, InvIS* mapping, int verbose){

    /* Get rank and size */
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* Open a file */
    char fullpath[256];
    snprintf(fullpath, sizeof(fullpath), "data/%s", filename);
    FILE* file = fopen(fullpath, "r");
    if(!file){
        if(verbose) printf("%sFile not found.\n", INV_MATRIX_VERBOSE);
        fclose(file);
        return NULL;
    } else {
        if(verbose) printf("%sReading a file \"%s\"...\n", INV_MATRIX_VERBOSE, fullpath);
    }

    /* Read the header */
    char* line = NULL;
    size_t line_len;
    int n, m, i, j;
    double val;
    getline(&line, &line_len, file);
    fscanf(file, "%i%i%i\n", &n, &n, &m);
    if(verbose) printf("%s%s", INV_MATRIX_VERBOSE, line+2);
    if(verbose) printf("%s%i rows, %i non-zeros.\n", INV_MATRIX_VERBOSE, n, m);

    /* Compute number of nonzeros under a reordering */
    int nnz[n];
    for(int k=0; k<n; k++) nnz[k] = 0;
    for(int k=0; k<m; k++){
        fscanf(file, "%i%i%lf\n", &i, &j, &val);
        nnz[mapping->index[i-1]]++;
        nnz[mapping->index[j-1]] += (i!=j);
    }

    /* Move the file pointer */
    fseek(file, 0, SEEK_SET);
    getline(&line, &line_len, file);
    fscanf(file, "%i%i%i\n", &n, &n, &m);

    /* Create a matrix */
    Mat Q;
    int Istart, Iend;
    int* nnz_local = nnz+mapping->offset[rank];
    MatCreate(comm, &Q);
    if(size>1){
        if(verbose) printf("%sCreating an MPI matrix on %i procs...\n", INV_MATRIX_VERBOSE, size);
        MatSetType(Q, MATMPIAIJ);
        MatSetSizes(Q, mapping->n_local[rank], mapping->n_local[rank], PETSC_DETERMINE, PETSC_DETERMINE);
        MatMPIAIJSetPreallocation(Q, 1, nnz_local, 1, nnz_local);
        // MatMPISBAIJSetPreallocation(Q, 1, 5, NULL, 2, NULL);    // must be changed if use MatMPIXXX
    } else {
        if(verbose) printf("%sCreating a SEQ matrix on %i proc...\n", INV_MATRIX_VERBOSE, size);
        MatSetType(Q, MATSEQSBAIJ);     // MATSEQAIJ?
        MatSetSizes(Q, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatSeqAIJSetPreallocation(Q, 1, nnz);
    }
    MatSetFromOptions(Q);
    MatSetUp(Q);
    MatGetOwnershipRange(Q, &Istart, &Iend);
    for(int k=0; k<m; k++){
        fscanf(file, "%i%i%lf\n", &i, &j, &val);
        if(mapping->part[i-1]==rank) MatSetValue(Q, mapping->index[i-1], mapping->index[j-1], val, INSERT_VALUES);
        if(mapping->part[j-1]==rank) MatSetValue(Q, mapping->index[j-1], mapping->index[i-1], val, INSERT_VALUES);
    }
    MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY);
    fclose(file);
    if(verbose) printf("%sMatrix assembled.\n", INV_MATRIX_VERBOSE);

    return Q;
}


Mat InvMatrixCreateLocalFromFile(MPI_Comm comm, char* filename, InvIS* mapping, int verbose){

    /* Get rank and size */
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* Open a file */
    char fullpath[256];
    snprintf(fullpath, sizeof(fullpath), "data/%s", filename);
    FILE* file = fopen(fullpath, "r");
    if(!file){
        if(verbose) printf("%sFile not found.\n", INV_MATRIX_VERBOSE);
        fclose(file);
        return NULL;
    } else {
        if(verbose) printf("%sReading a file \"%s\"...\n", INV_MATRIX_VERBOSE, fullpath);
    }

    /* Read the header */
    char* line = NULL;
    size_t line_len;
    int n, m, i, j;
    double val;
    getline(&line, &line_len, file);
    fscanf(file, "%i%i%i\n", &n, &n, &m);
    if(verbose) printf("%s%s", INV_MATRIX_VERBOSE, line+2);
    if(verbose) printf("%s%i rows, %i non-zeros.\n", INV_MATRIX_VERBOSE, n, m);

    /* Compute number of nonzeros under a reordering */
    int nnz[n];
    for(int k=0; k<n; k++) nnz[k] = 0;
    for(int k=0; k<m; k++){
        fscanf(file, "%i%i%lf\n", &i, &j, &val);
        if(mapping->part[i-1]==rank && mapping->part[j-1]==rank){
            nnz[mapping->index[i-1]]++;
            nnz[mapping->index[j-1]] += (i!=j);
        }
    }

    /* Move the file pointer */
    fseek(file, 0, SEEK_SET);
    getline(&line, &line_len, file);
    fscanf(file, "%i%i%i\n", &n, &n, &m);

    /* Create a matrix */
    if(verbose) printf("%sCreating a SEQ matrix on %i proc...\n", INV_MATRIX_VERBOSE, size);
    Mat Q;
    int Istart, Iend;
    MatCreate(PETSC_COMM_SELF, &Q);
    MatSetType(Q, MATSEQSBAIJ);     // MATSEQAIJ?
    MatSetSizes(Q, PETSC_DECIDE, PETSC_DECIDE, mapping->n_local[rank], mapping->n_local[rank]);
    MatSeqAIJSetPreallocation(Q, 1, nnz);
    MatSetFromOptions(Q);
    MatSetUp(Q);
    MatGetOwnershipRange(Q, &Istart, &Iend);
    for(int k=0; k<m; k++){
        fscanf(file, "%i%i%lf\n", &i, &j, &val);
        if(mapping->part[i-1]==rank && mapping->part[j-1]==rank){
            MatSetValue(Q, mapping->index[i-1], mapping->index[j-1], val, INSERT_VALUES);
            MatSetValue(Q, mapping->index[j-1], mapping->index[i-1], val, INSERT_VALUES);
        }
    }
    MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY);
    fclose(file);
    if(verbose) printf("%sMatrix assembled.\n", INV_MATRIX_VERBOSE);

    return Q;
}


Mat InvMatrixFactorLocal(Mat Q, int verbose){

    if(verbose) printf("%sFactoring the local matrix...\n", INV_MATRIX_VERBOSE);
    Mat L;
    IS isrow, iscol;
    MatSetOption(Q, MAT_SYMMETRIC, PETSC_TRUE);
    MatGetOrdering(Q, MATORDERINGNATURAL, &isrow, &iscol);
    // MatDuplicate(Q, MAT_SHARE_NONZERO_PATTERN, &L);
    MatGetFactor(Q, MATSOLVERMUMPS, MAT_FACTOR_CHOLESKY, &L);
    MatCholeskyFactorSymbolic(L, Q, isrow, NULL);
    MatCholeskyFactorNumeric(L, Q, NULL);
    ISDestroy(&isrow);
    ISDestroy(&iscol);
    if(verbose) printf("%sFactoring completed.\n", INV_MATRIX_VERBOSE);

    return L;
}

Mat InvMatrixCreateLocalSubmatrix(MPI_Comm comm, Mat Q, InvIS* mapping, InvIS* mapping_sub, int verbose){

    /* Get rank and size */
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* Get local indices */
    if(verbose) printf("%sCreating a submatrix...\n", INV_MATRIX_VERBOSE);
    int temp[mapping_sub->n_local[rank]], i_local = 0, j_local = 0;
    for(int i=0; i<mapping_sub->n_local[rank]; i++) temp[i] = 0;
    for(int i=0; i<mapping->n_vert; i++){
        if(mapping->part[i]==rank && mapping_sub->part[i]==rank){
            temp[i_local] = j_local;
            i_local++;
        }
        j_local += mapping->part[i]==rank;
    }

    /* Get local submatrix */
    Mat Q_sub;
    IS isrow;
    ISCreateGeneral(PETSC_COMM_SELF, mapping_sub->n_local[rank], temp, PETSC_COPY_VALUES, &isrow);
    MatCreateSubMatrix(Q, isrow, isrow, MAT_INITIAL_MATRIX, &Q_sub);
    if(verbose) printf("%sSubmatrix created.\n", INV_MATRIX_VERBOSE);

    return Q_sub;
}


Mat InvMatrixCreateRW2(MPI_Comm comm, int n, double alpha, double tau, int cyclic, int verbose){

    int size;
    MPI_Comm_size(comm, &size);

    Mat Q;
    int Istart, Iend;
    MatCreate(comm, &Q);
    if(size>1){
        if(verbose) printf("%sCreating an MPI matrix on %i procs...\n", INV_MATRIX_VERBOSE, size);
        MatSetType(Q, MATMPIAIJ);
        MatSetSizes(Q, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatMPISBAIJSetPreallocation(Q, 1, 5, NULL, 2, NULL);
    } else {
        if(verbose) printf("%sCreating a SEQ matrix on %i proc...\n", INV_MATRIX_VERBOSE, size);
        MatSetType(Q, MATSEQSBAIJ);     // MATSEQAIJ?
        MatSetSizes(Q, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatSeqSBAIJSetPreallocation(Q, 1, 5, NULL);
    }
    MatSetFromOptions(Q);
    MatSetUp(Q);
    MatGetOwnershipRange(Q, &Istart, &Iend);
    for(int i=Istart; i<Iend; i++){
        if(cyclic){
            MatSetValue(Q, i, (i-2+n)%n, 1 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, (i-1+n)%n, -4 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, (i+0+n)%n, 6 * alpha + tau, INSERT_VALUES);
            MatSetValue(Q, i, (i+1+n)%n, -4 * alpha, INSERT_VALUES);
            MatSetValue(Q, i, (i+2+n)%n, 1 * alpha, INSERT_VALUES);
        }
        else if(i==0){
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

