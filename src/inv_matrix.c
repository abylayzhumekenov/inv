#include <petsc.h>
#include <petscblaslapack.h>

#include "inv_matrix.h"
#include "inv_pardiso.h"


PetscErrorCode MatCreateLoad(MPI_Comm comm, MatType type, PetscInt m, PetscInt n, PetscInt M, PetscInt N, const char filename[], Mat* A){

    int petsc_decide = (m==PETSC_DECIDE && n==PETSC_DECIDE && M==PETSC_DETERMINE && N==PETSC_DETERMINE);
    PetscViewer viewer;
    MatCreate(comm, A);
    MatSetType(*A, type);
    if(!petsc_decide) MatSetSizes(*A, m, n, M, N);
    PetscViewerBinaryOpen(comm, filename, FILE_MODE_READ, &viewer);
    MatLoad(*A, viewer);

    PetscViewerDestroy(&viewer);
    return 0;
}


PetscErrorCode MatGetSizeLoad(const char filename[], int* m, int* n){

    FILE* file = fopen(filename, "rb");
    unsigned int header[4];
    fread(header, 4, 4, file);
    *m = htobe32(header[1]);
    *n = htobe32(header[2]);
    fclose(file);

    return 0;
}


PetscErrorCode MatCreateSubMatrixSeq(Mat A, IS isrow, IS iscol, Mat* B){

    Mat *BB[1];
    MatCreateSubMatrices(A, 1, &isrow, &iscol, MAT_INITIAL_MATRIX, BB);
    *B = *BB[0];
    
    return 0;
}


PetscErrorCode MatMPIAIJKron(Mat A, Mat B, Mat* C){

    int nt_local, mt_local, ns, ms;
    Mat *A_seq[1], C_seq;
    IS isrow, iscol;

    MatGetOwnershipIS(A, &isrow, &iscol);
    MatCreateSubMatrices(A, 1, &isrow, &iscol, MAT_INITIAL_MATRIX, A_seq);
    MatSeqAIJKron(*A_seq[0], B, MAT_INITIAL_MATRIX, &C_seq);
    MatGetSize(*A_seq[0], &nt_local, &mt_local);
    MatGetSize(B, &ns, &ms);
    MatCreateMPIMatConcatenateSeqMat(PETSC_COMM_WORLD, C_seq, nt_local*ms, MAT_INITIAL_MATRIX, C);

    MatDestroyMatrices(1, A_seq);
    MatDestroy(&C_seq);
    return 0;
}


PetscErrorCode MatDenseInvertLapack(Mat A){

    int rank, size, n;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MatGetLocalSize(A, &n, &n);

    double *a_array, work[n];
    int i_array[n], info;
    MatDenseGetArray(A, &a_array);
    if(rank==size-1){
        dgetrf_(&n, &n, a_array, &n, i_array, &info);
        dgetri_(&n, a_array, &n, i_array, work, &n, &info);
    }
    MatDenseRestoreArray(A, &a_array);

    return 0;
}


PetscErrorCode MatDensePointwiseMult(Mat A, Mat B){

    double *a_array, *b_array;
    int n, m;
    MatGetLocalSize(A, &n, NULL);
    MatGetSize(A, NULL, &m);
    MatDenseGetArray(A, &a_array);
    MatDenseGetArray(B, &b_array);
    for(int i=0; i<n*m; i++) a_array[i] *= b_array[i];
    MatDenseRestoreArray(A, &a_array);
    MatDenseRestoreArray(B, &b_array);

    return 0;
}

PetscErrorCode MatSeqAIJInvertMUMPS(Mat A, Mat* B, int verbose){

    /* Declare the variables */
    if(verbose) printf("\tDirect solver:\t\tMUMPS\n");
    int n;
    Mat L;
    IS is_factor, is_factor_col;

    /* Factor the matrix */
    if(verbose) printf("\tFactoring the matrix...\n");
    MatGetSize(A, &n, &n);
    MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);
    MatGetOrdering(A, MATORDERINGNATURAL, &is_factor, &is_factor_col);
    MatGetFactor(A, MATSOLVERMUMPS, MAT_FACTOR_CHOLESKY, &L);
    MatCholeskyFactorSymbolic(L, A, is_factor, NULL);
    MatCholeskyFactorNumeric(L, A, NULL);

    /* Sparse inversion */
    if(verbose) printf("\tSolving sparse rhs...\n");
    MatCreate(PETSC_COMM_SELF, B);
    MatSetSizes(*B, n, n, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetType(*B, MATSEQAIJ);
    MatSetFromOptions(*B);
    MatSetUp(*B);
    for(int i=0; i<n; i++) MatSetValue(*B, i, i, 1.0, INSERT_VALUES);
    MatAssemblyBegin(*B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*B, MAT_FINAL_ASSEMBLY);
    MatMumpsGetInverseTranspose(L, *B);

    /* Clean up */
    if(verbose) printf("\tClean up...\n");
    MatDestroy(&L);
    ISDestroy(&is_factor);
    ISDestroy(&is_factor_col);

    return 0;
}

PetscErrorCode MatSeqAIJInvertPARDISO(Mat A, Mat* B, int verbose){

        /* Declare variables */
        if(verbose) printf("\tDirect solver:\t\tPARDISO\n");
        int n_temp, n, nnz;
        const int *ia_temp, *ja_temp;
        int *ia, *ja;
        double *val;
        PetscBool done_temp;

        /* Copy the matrix data */
        if(verbose) printf("\tCopying the matrix data...\n");
        MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);
        MatConvert(A, MATSBAIJ, MAT_INITIAL_MATRIX, B);
        MatSeqAIJGetArray(*B, &val);
        MatGetRowIJ(*B, 0, PETSC_FALSE, PETSC_FALSE, &n_temp, &ia_temp, &ja_temp, &done_temp);
        n = n_temp;
        nnz = ia_temp[n];
        ia = calloc(sizeof(int), n+1);
        ja = calloc(sizeof(int), nnz);
        memcpy(ia, ia_temp, sizeof(int)*(n+1));
        memcpy(ja, ja_temp, sizeof(int)*nnz);
        MatRestoreRowIJ(*B, 0, PETSC_FALSE, PETSC_FALSE, &n_temp, &ia_temp, &ja_temp, &done_temp);

        /* PARDISO control parameters */
        int      nrhs = 1;  /* Number of rhs */
        int      mtype = -2;/* Symmetric matrix */
        void    *pt[64];    /* Internal pointer */
        int      iparm[64];
        double   dparm[64];
        int      maxfct, mnum, phase, error = 0, msglvl, solver = 0;
        int      num_procs = 1; /* Default number of threads */
        char    *var;       /* Enironment variable */
        double   ddum;      /* Double dummy. */
        int      idum;      /* Integer dummy. */

        /* Initialize PARDISO */
        if(verbose) printf("\tInitializing PARDISO...\n");
        omp_set_num_threads(1);
        var = getenv("OMP_NUM_THREADS");
        if(var != NULL) sscanf(var, "%d", &num_procs);
        pardisoinit (pt, &mtype, &solver, iparm, dparm, &error);
        iparm[2]  = num_procs;
        maxfct = 1;		    /* Maximum number of numerical factorizations */
        mnum   = 1;         /* Which factorization to use */
        msglvl = 0;         /* Print statistical information */
        error  = 0;         /* Initialize error flag */

        /* Convert to 1-based Fortran indexing */
        for(int i=0; i<n+1; i++) ia[i] += 1;
        for(int i=0; i<nnz; i++) ja[i] += 1;

        /* Reordering and symbolic factorization.  */
        /* This step also allocates all memory that is necessary for the factorization. */
        if(verbose) printf("\tSymbolic factorization...\n");
        phase = 11;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
                &n, val, ia, ja, &idum, &nrhs, iparm, 
                &msglvl, &ddum, &ddum, &error, dparm);

        /* Numerical factorization */
        if(verbose) printf("\tNumerical factorization...\n");
        phase = 22;
        iparm[32] = 0;  /* compute determinant */
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
                &n, val, ia, ja, &idum, &nrhs, iparm, 
                &msglvl, &ddum, &ddum, &error, dparm);

        /* Selected inversion */
        if(verbose) printf("\tSelected inversion...\n");
        phase = -22;
        iparm[35] = 1;  /* do not overwrite internal factor L */
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
                &n, val, ia, ja, &idum, &nrhs, iparm, 
                &msglvl, NULL, NULL, &error, dparm);

        // /* Convert back to 0-based C indexing (for freeing?) */
        // for(int i=0; i<n+1; i++) ia[i] -= 1;
        // for(int i=0; i<nnz; i++) ja[i] -= 1;

        /* Finalize PARDISO */
        if(verbose) printf("\tClean up...\n");
        phase = -1;     /* Release internal memory. */
        pardiso(pt, &maxfct, &mnum, &mtype, &phase,
                &n, &ddum, ia, ja, &idum, &nrhs, iparm,
                &msglvl, &ddum, &ddum, &error, dparm);
        
        /* Free space and restore array */
        free(ia);
        free(ja);
        MatSeqAIJRestoreArray(*B, &val);
        
        return 0;
}