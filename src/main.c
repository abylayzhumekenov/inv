#include <stdlib.h>
#include <stdio.h>

#include <pcg_variants.h>
#include <petsc.h>
#include <metis.h>

#include "inv_sampler.h"
#include "inv_random.h"
#include "inv_input.h"
#include "inv_matrix.h"
#include "inv_graph.h"
#include "inv_is.h"
#include "inv_solver.h"


int main(int argc, char** argv){

    /* Set parameters from the options */
    int n, max_niter = 1000, n_samples = 1000, n_neighbors = 2, verbose = 0, verbose_sampler = 0;
    for(int i=0; i<argc; i++){
        if(!strcmp(argv[i], "-ns")){
            n_samples = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nmax")){
            max_niter = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-nn")){
            n_neighbors = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-v")){
            verbose = atoi(argv[i+1]);
        } else if(!strcmp(argv[i], "-vs")){
            verbose_sampler = atoi(argv[i+1]);
        }
        argv[i][0] = 0;
    }

    // ----------------------------------------------------------------------

    /* Declare variables */
    Vec *x, *z, d_sample, d_base, v;
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
    InputCSR* input = input_read("Qmm", verbose);

    /* Create a sequential graph */
    // Graph* graph = graph_create_rw2(100, verbose);
    Graph* graph = graph_create_from_input(input, verbose);

    /* Separate the graph recursively */
    /* Or perform k-way partitioning */
    /* Then no need to treat separator differently, and all ranks are balanced. */
    // graph = graph_separate(graph, verbose);
    graph = graph_partition(graph, world_size, verbose);
    is_part = is_create_from_partition(PETSC_COMM_SELF, graph, world_rank, verbose);

    // ----------------------------------------------------------------------

    /* Create a sequential matrix */
    // Q = matrix_create_rw2(PETSC_COMM_SELF, 100, 1.0, 1e-4, verbose);
    Q = matrix_create_from_input(PETSC_COMM_SELF, input, verbose);
    MatGetSize(Q, &n, &n);
    ISCreateStride(PETSC_COMM_SELF, n, 0, 1, &is_x);

    /* Sample */
    /* For now, only a sequential sampler works. */
    /* Need a parallel symmetric preconditioner for parallel sampling. */
    ksp_sampler = sampler_ksp(PETSC_COMM_SELF, Q, max_niter, verbose);
    z = sampler_std_normal(PETSC_COMM_SELF, n, n_samples, verbose);
    x = sampler_gmrf(ksp_sampler, z, n_samples, verbose, verbose_sampler);

    // ----------------------------------------------------------------------

    /* Create a submatrix Q[a,a] */
    graph = graph_neighborhood(graph, world_rank, n_neighbors, verbose);
    is_a = is_create_from_partition(PETSC_COMM_SELF, graph, world_rank, verbose);
    Q_a = matrix_create_submatrix(Q, is_a, verbose);

    /* Grow once more to include enclosing separator and create a submatrix Q[as,as] */
    graph = graph_neighborhood(graph, world_rank, 1, verbose);
    is_as = is_create_from_partition(PETSC_COMM_SELF, graph, world_rank, verbose);
    // Q_as = matrix_create_submatrix(Q, is_as, verbose);

    // ----------------------------------------------------------------------

    /* Compute the sample contribution */
    ISDifference(is_as, is_a, &is_s);
    ISDifference(is_x, is_s, &is_x);
    ksp_solver = solver_ksp(PETSC_COMM_SELF, Q_a, max_niter, verbose);
    d_sample = solver_sample_contribution(ksp_solver, Q, x, is_x, is_a, n_samples, verbose);

    /* Solve the base case Q[a,a]^-1 and add */
    d_base = solver_base_case(PETSC_COMM_SELF, Q_a, verbose);
    VecAXPY(d_sample, 1.0, d_base);

    /* Assemble the solution */
    v = solver_assemble_solution(d_sample, is_part, is_a, n, verbose);

    // ----------------------------------------------------------------------

    /* Output */
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "data/x", FILE_MODE_WRITE, &viewer);
    VecView(v, viewer);

    // ----------------------------------------------------------------------

    /* Free objects */
    input_destroy(input);
    graph_destroy(graph);

    /* Free PETSc objects */
    for(int i=0; i<n_samples; i++){
        VecDestroy(&z[i]);
        VecDestroy(&x[i]);
    }
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
