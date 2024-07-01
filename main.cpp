#include "laplace_solver.hpp"
#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "muparser_fun.hpp"


int main(int argc, char* argv[]) {

    Timings::Chrono     clock1;
    clock1.start();

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::ifstream f("parameters.json");
    json paras = json::parse(f);

    std::string funString = paras.value("fun","");
    MuparserFun Fun(funString);

    bool from_script = false;

    // Verifica se è stato passato l'argomento "test"
    
    int n = paras.value("n", 10);
    int num_threads = paras.value("num_threads", 4);


    if (argc > 2) {
    n = atoi(argv[1]);
    num_threads = atoi(argv[2]);
        if (argc == 4){
        if (std::string(argv[3]) == "test") {
            from_script = true;
        }}
    }

    omp_set_num_threads(num_threads);

    double h = 1.0 / (n - 1);
    double tol = paras.value("tol", 1e-5);
    int max_iter = paras.value("max_it",1000000);

    double* U = (double*)calloc(n * n, sizeof(double));
    double* U_new = (double*)calloc(n * n, sizeof(double));

    

    laplace_solver(U, U_new, n, h, f_function, rank, size, tol, max_iter, num_threads, from_script);


    
    free(U);
    free(U_new);

    MPI_Finalize();
    
    clock1.stop();
    
    if(rank == 0)
    std::cout << clock1;
    

    return 0;
}

