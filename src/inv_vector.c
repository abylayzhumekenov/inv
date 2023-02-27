#include "inv_vector.h"


PetscErrorCode VecCreateLoad(MPI_Comm comm, VecType type, PetscInt n, PetscInt N, const char filename[], Vec* v){

    int petsc_decide = (n==PETSC_DECIDE && N==PETSC_DETERMINE);
    PetscViewer viewer;
    VecCreate(comm, v);
    VecSetType(*v, type);
    if(!petsc_decide) VecSetSizes(*v, n, N);
    PetscViewerBinaryOpen(comm, filename, FILE_MODE_READ, &viewer);
    VecLoad(*v, viewer);

    PetscViewerDestroy(&viewer);
    return 0;
}


PetscErrorCode VecCreateSubVector(Vec x, IS is, Vec* y){

    int n;
    VecScatter scatter;
    ISGetLocalSize(is, &n);
    VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, y);
    VecScatterCreate(x, is, *y, NULL, &scatter);
    VecScatterBegin(scatter, x, *y, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter, x, *y, INSERT_VALUES, SCATTER_FORWARD);

    VecScatterDestroy(&scatter);
    return 0;
}

PetscErrorCode VecCreateSubVectorSeq(Vec x, IS is, Vec* y){

    int n;
    Vec v;
    double* v_array;
    ISGetLocalSize(is, &n);
    VecCreateSubVector(x, is, &v);
    VecGetArray(v, &v_array);
    VecCreateSeq(PETSC_COMM_SELF, n, y);
    for(int i=0; i<n; i++) VecSetValue(*y, i, v_array[i], INSERT_VALUES);
    
    VecDestroy(&v);
    return 0;
}


PetscErrorCode VecCreateSuperVector(Vec x, Vec* y){

    int n;
    Vec y_local;
    VecGetLocalSize(x, &n);
    VecCreateMPI(PETSC_COMM_WORLD, n, PETSC_DETERMINE, y);
    VecCreateLocalVector(*y, &y_local);
    VecGetLocalVector(*y, y_local);
    VecCopy(x, y_local);
    VecRestoreLocalVector(*y, y_local);

    return 0;
}