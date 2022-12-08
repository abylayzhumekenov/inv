#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <petsc.h>
#include <mpi.h>

#include "inv_is.h"
#include "inv_graph.h"
#include "inv_matrix.h"
#include "inv_sampler.h"
#include "inv_solver.h"


int main(int argc, char **argv){

    /* Set parameters from the options */
    char *filein = "Qmm", fileout[256];
    int max_niter = 1000, n_sample = 1000, n_neighbor = 1, verbose = 0, verbose_s = 0, n_part = 1;
    for(int i=0; i<argc; i++){
        if(!strcmp(argv[i], "-fin")){
            char filetemp[128];
            strcpy(filetemp, argv[i+1]);
            filein = filetemp;
            strcpy(fileout, "data/x_");
            strcat(fileout, argv[i+1]);
        } else if(!strcmp(argv[i], "-ns")){
            n_sample = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nmax")){
            max_niter = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nn")){
            n_neighbor = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-v")){
            verbose = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-vs")){
            verbose_s = atoi(argv[i+1]);
        }
        argv[i] = 0;
    }


    /* Initialize MPI */
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    verbose = (!rank) && verbose, n_part = size;


    /* Initialize PETSc */
    PetscInitialize(&argc, &argv, (char*)0, (char*)0);
    MPI_Comm comm = PETSC_COMM_WORLD;
    PetscViewer viewer = PETSC_VIEWER_STDOUT_WORLD;


    /* Read a graph */
    InvGraph* graph = InvGraphCreateFromFile(filein, verbose);
    graph = InvGraphPartition(graph, n_part, verbose);


    /* Create a mappings */
    InvIS* itog = InvISCreateInitialToGlobal(graph, verbose);
    InvIS* gtoi = InvISCreateGlobalToInitial(graph, verbose);


    /* Read a matrix */
    Mat Q = InvMatrixCreateFromFile(comm, filein, itog, verbose);


    /* Compute halo regions */
    graph = InvGraphNeighborhood(graph, rank, n_neighbor, verbose);
    InvIS* itol = InvISCreateInitialToLocal(graph, verbose);
    InvIS* ltoi = InvISCreateLocalToInitial(graph, rank, verbose);
    graph = InvGraphNeighborhood(graph, rank, 1, verbose);
    InvIS* itol2 = InvISCreateInitialToLocal(graph, verbose);
    InvIS* ltoi2 = InvISCreateLocalToInitial(graph, rank, verbose);


    /* Read ghost matrices */
    Mat Q_local = InvMatrixCreateLocalFromFile(comm, filein, itol2, verbose);
    Mat Q_sub = InvMatrixCreateLocalSubmatrix(comm, Q_local, itol2, itol, verbose);


    /* Solve the base case */
    Mat L = InvMatrixFactorLocal(Q_sub, verbose);
    Vec d = InvSolverBaseCase(L, verbose);


    /* Create a sampler KSP */
    pcg64_random_t rng = InvRandomCreate(0, 0, 0);
    KSP ksp_sampler = InvSamplerCreateKSP(comm, Q, max_niter, verbose);
    Vec z = NULL, x = NULL, y, w, v;


    double log;
    /* Combine samples and direct solution */
    for(int i=0; i<n_sample; i++){
        if(verbose && verbose_s) printf("SAMPLER:\tSample #%i\n", i);
        z = InvSamplerStdNormal(comm, z, &rng, itog, verbose && verbose_s);
        x = InvSamplerGMRF(ksp_sampler, x, z, verbose && verbose_s);
        y = InvSolverMultQx(comm, Q_local, x, itol2, itol, itog, verbose && verbose_s);
        w = InvSolverSampleSq(L, y, verbose && verbose_s);
        VecAXPY(d, 1.0/n_sample, w);

        VecDestroy(&y);
        VecDestroy(&w);
        VecDestroy(&z);
        VecDestroy(&x);
        PetscMallocGetCurrentUsage(&log);
        if(!rank) printf("----------------- %f \n", log);
        // free z, x, y, w?
        // or rewrite z, x, y, w...
    }


    /* Assemble the solution */
    v = InvSolverAssembleSolution(comm, d, itol, itog, gtoi, verbose);


    /* Write the solution */
    PetscViewerBinaryOpen(comm, fileout, FILE_MODE_WRITE, &viewer);
    VecView(v, viewer);


    /* Clean up */
    VecDestroy(&z);
    VecDestroy(&x);
    VecDestroy(&y);
    VecDestroy(&w);
    VecDestroy(&v);
    MatDestroy(&Q);
    MatDestroy(&Q_local);
    MatDestroy(&Q_sub);
    MatDestroy(&L);
    KSPDestroy(&ksp_sampler);
    PetscViewerDestroy(&viewer);
    InvGraphDestroy(graph, verbose);
    InvISDestroy(itog, verbose);
    InvISDestroy(gtoi, verbose);
    InvISDestroy(itol, verbose);
    InvISDestroy(ltoi, verbose);
    InvISDestroy(itol2, verbose);
    InvISDestroy(ltoi2, verbose);


    /* Finalize */
    PetscFinalize();
    MPI_Finalize();

    return 0;
}