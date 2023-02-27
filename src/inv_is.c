#include "inv_is.h"


PetscErrorCode ISCreateLoad(MPI_Comm comm, const char filename[], IS* is){

    PetscViewer viewer;
    ISCreate(comm, is);
    PetscViewerBinaryOpen(comm, filename, FILE_MODE_READ, &viewer);
    ISLoad(*is, viewer);

    PetscViewerDestroy(&viewer);
    return 0;
}


PetscErrorCode MatGetOwnershipISRows(MPI_Comm comm, Mat A, IS* is){

    int istart, iend, n;
    MatGetOwnershipRange(A, &istart, &iend);
    MatGetLocalSize(A, &n, NULL);
    ISCreateStride(comm, n, istart, 1, is);
    ISToGeneral(*is);

    return 0;
}


PetscErrorCode ISCreateBlockIS(MPI_Comm comm, IS is, int bs, IS* is_block){

    int n;
    const int *is_array;
    ISGetLocalSize(is, &n);
    ISGetIndices(is, &is_array);
    ISCreateBlock(comm, bs, n, is_array, PETSC_COPY_VALUES, is_block);
    ISRestoreIndices(is, &is_array);
    
    return 0;
}