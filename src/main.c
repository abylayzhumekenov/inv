#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <pcg_variants.h>
#include <petsc.h>
#include <metis.h>
#include <sbaij.h>

#include "inv_sampler.h"
#include "inv_random.h"
#include "inv_input.h"
#include "inv_matrix.h"
#include "inv_graph.h"
#include "inv_is.h"
#include "inv_solver.h"


int main(int argc, char** argv){

    /* Set parameters from the options */
    char *filein = "Qmm", fileout[256];
    int max_niter = 1000, n_samples = 1000, n_neighbors = 2, verbose = 0;
    for(int i=0; i<argc; i++){
        if(!strcmp(argv[i], "-fin")){
            char filetemp[128];
            strcpy(filetemp, argv[i+1]);
            filein = filetemp;
            strcpy(fileout, "data/x_");
            strcat(fileout, argv[i+1]);
        } else if(!strcmp(argv[i], "-ns")){
            n_samples = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nmax")){
            max_niter = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nn")){
            n_neighbors = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-v")){
            verbose = atoi(argv[i+1]);
        }
        argv[i] = 0;
    }
    

    // ----------------------------------------------------------------------

    /* Declare variables */
    int n, n_a;
    Vec x, z, d, d_sample, d_base, v;
    Mat Q, Q_a;
    KSP ksp_sampler, ksp_solver;
    IS is_part, is_x, is_a, is_as, is_s;
    PetscViewer viewer;

    /* Initialize MPI */
    MPI_Init(NULL, NULL);

    /* Create communicators */
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    verbose = verbose && !world_rank;

    /* Initialize PETSc */
    PETSC_COMM_WORLD = MPI_COMM_WORLD;
    PetscInitialize(&argc, &argv, (char*)0, (char*)0);

    // ----------------------------------------------------------------------
    
    /* Read the input file */
    InputCSR* input = input_read(filein, verbose);

    /* Create a sequential graph */
    // Graph* graph = graph_create_rw2(10, verbose);
    Graph* graph = graph_create_from_input(input, verbose);

    /* Separate the graph recursively OR perform k-way partitioning */
    /* Then no need to treat separator differently, and all ranks are balanced. */
    // graph = graph_separate(graph, verbose);
    graph = graph_partition(graph, world_size, verbose);
    is_part = is_create_from_partition(PETSC_COMM_SELF, graph, world_rank, verbose);

    // ----------------------------------------------------------------------

    /* Create a sequential matrix */
    // Q = matrix_create_rw2(PETSC_COMM_SELF, 10, 1.0, 1e-4, verbose);
    Q = matrix_create_from_input(PETSC_COMM_SELF, input, verbose);
    MatGetSize(Q, &n, &n);
    ISCreateStride(PETSC_COMM_SELF, n, 0, 1, &is_x);

    /* Create a submatrix Q[a,a] */
    graph = graph_neighborhood(graph, world_rank, n_neighbors, verbose);
    is_a = is_create_from_partition(PETSC_COMM_SELF, graph, world_rank, verbose);
    Q_a = matrix_create_submatrix(Q, is_a, verbose);

    /* Grow once more to include enclosing separator and create a submatrix Q[as,as] */
    graph = graph_neighborhood(graph, world_rank, 1, verbose);
    is_as = is_create_from_partition(PETSC_COMM_SELF, graph, world_rank, verbose);
    ISDifference(is_as, is_a, &is_s);
    ISDifference(is_x, is_s, &is_x);

    // ----------------------------------------------------------------------

    /* Prepare KSP solvers and vectors */
    ksp_sampler = sampler_ksp(PETSC_COMM_SELF, Q, max_niter, verbose);
    ksp_solver = solver_ksp(PETSC_COMM_SELF, Q_a, max_niter, verbose);
    ISGetSize(is_a, &n_a);
    VecCreate(PETSC_COMM_SELF, &d);
    VecSetType(d, VECSEQ);
    VecSetSizes(d, PETSC_DECIDE, n_a);
    VecSetFromOptions(d);
    VecSetUp(d);
    VecSet(d, 0.0);

    /* Create a random number generator */
    pcg64_random_t rng;
    pcg64_srandom_r(&rng, time(NULL), 0);
    // pcg64_srandom_r(&rng, time(NULL), (intptr_t)&rng);  // completely unpredictable seeds

    /* Compute the variance term */
    for(int k=0; k<n_samples; k++){

        /* Sample */
        /* For now, only a sequential sampler works. */
        /* Need a parallel symmetric preconditioner for parallel sampling. */
        z = sampler_std_normal(PETSC_COMM_SELF, &rng, n, verbose, k);
        x = sampler_gmrf(ksp_sampler, z, verbose, k);

        /* Compute the sample contribution */
        d_sample = solver_sample_contribution(ksp_solver, Q, x, is_x, is_a, verbose);
        VecAXPY(d, 1.0 / n_samples, d_sample);
    }

    /* Solve the base case Q[a,a]^-1 and add */
    d_base = solver_base_case(PETSC_COMM_SELF, Q_a, verbose);
    VecAXPY(d, 1.0, d_base);

    /* Assemble the solution */
    v = solver_assemble_solution(d, is_part, is_a, n, verbose);

    // ----------------------------------------------------------------------

    /* Output */
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileout, FILE_MODE_WRITE, &viewer);
    VecView(v, viewer);

    // ----------------------------------------------------------------------

    /* Free objects */
    input_destroy(input);
    graph_destroy(graph);

    /* Free PETSc objects */
    VecDestroy(&z);
    VecDestroy(&x);
    VecDestroy(&d);
    VecDestroy(&d_sample);
    VecDestroy(&d_base);
    VecDestroy(&v);
    MatDestroy(&Q);
    MatDestroy(&Q_a);
    KSPDestroy(&ksp_sampler);
    KSPDestroy(&ksp_solver);
    ISDestroy(&is_part);
    ISDestroy(&is_x);
    ISDestroy(&is_a);
    ISDestroy(&is_as);
    ISDestroy(&is_s);
    PetscViewerDestroy(&viewer);

    // ----------------------------------------------------------------------    

    /* Finalize PETSc*/
    PetscFinalize();

    /* Finalize MPI */
    MPI_Finalize();

    return 0;
}
