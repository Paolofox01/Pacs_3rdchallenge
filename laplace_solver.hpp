#ifndef LAPLACE_SOLVER_H
#define LAPLACE_SOLVER_H

#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <chrono.hpp>
#include "json.hpp"
using json = nlohmann::json;
#include "muparser_fun.hpp"


void laplace_solver(double* U, double* U_new, int n, double h, double (*f)(double, double), int rank, int size, double tol, int max_iter, int num_threads, bool from_script);
void initialize_grid(double* U, int n, int rank, int size, int start_row, int end_row);
void update_grid(double* U, double* U_new, int n, double h, double (*f)(double, double), int rank, int size, int start_row, int end_row);
bool check_convergence(double* U, double* U_new, int n, double tol, int rank, int size, double h, int start_row, int end_row);
double f_function(double x, double y);
void export_to_vtk(double* U, int n, int rank, int size, int num_threads, bool from_script);




#endif
