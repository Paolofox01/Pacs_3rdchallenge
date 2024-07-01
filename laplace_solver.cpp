#include "laplace_solver.hpp" //@note it is called laplace_solver.hpp!


/**
 * @brief Risolve l'equazione di Laplace su una griglia quadrata.
 * 
 * @param U Puntatore all'array che contiene i valori della griglia
 * @param U_new Puntatore all'array per i nuovi valori della griglia
 * @param n Dimensione della griglia
 * @param h Passo della griglia
 * @param f Funzione di forzamento
 * @param rank Rank del processo MPI
 * @param size Numero totale di processi MPI
 * @param tol Tolleranza per la convergenza
 * @param max_iter Numero massimo di iterazioni
 * @param num_threads Numero di thread OpenMP
 * @param from_script Se il codice è stato lanciato da uno script
 */
//@note why pointers. In C++ we can use references, which are more readable and safer
// also use std::function instead a function pointer 
// Why not making a class, with a method solve, and the data as members?
void laplace_solver(double* U, double* U_new, int n, double h, double (*f)(double, double), int rank, int size, double tol, int max_iter, int num_threads, bool from_script){
    
    int rem = n % size;
    int local_rows = (rank < rem) ? (n / size) + 1 :  n / size;
    
    int start_row = (rank < rem) ? rank * (local_rows) :  rem * (local_rows + 1) + (rank - rem) * local_rows;
    int end_row = start_row + local_rows - 1;
    
    
    initialize_grid(U, n, rank, size, start_row, end_row);

    bool converged = false;
    for (int iter = 0; iter < max_iter; iter++) {
        update_grid(U, U_new, n, h, f_function, rank, size, start_row, end_row);
        converged = check_convergence(U, U_new, n, tol, rank, size, h, start_row, end_row);
        if (converged) break;
        std::swap(U, U_new);
    }

    MPI_Gather(U + start_row * n, local_rows * n, MPI_DOUBLE, U + start_row * n, local_rows * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (rank == 0){ 
        export_to_vtk(U, n, rank, size, num_threads, from_script);}
    


}

/**
 * @brief Inizializza la griglia con condizioni iniziali e di bordo.
 * 
 * @param U Puntatore all'array che contiene i valori della griglia
 * @param n Dimensione della griglia
 * @param rank Rank del processo MPI
 * @param size Numero totale di processi MPI
 * @param start_row Indice di inizio delle righe locali
 * @param end_row Indice di fine delle righe locali
 */
void initialize_grid(double* U, int n, int rank, int size, int start_row, int end_row) {

    #pragma omp parallel for
    for (int i = start_row; i < end_row; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == 0 || i == n-1 || j == 0 || j == n-1) {
                U[i * n + j] = 0.0; // boundary condition
            } else {
                U[i * n + j] = 0.0; // initial guess
            }
        }
    }
}

/**
 * @brief Aggiorna la griglia con il metodo di Jacobi.
 * 
 * @param U Puntatore all'array che contiene i valori della griglia
 * @param U_new Puntatore all'array per i nuovi valori della griglia
 * @param n Dimensione della griglia
 * @param h Passo della griglia
 * @param f Funzione di forzamento
 * @param rank Rank del processo MPI
 * @param size Numero totale di processi MPI
 * @param start_row Indice di inizio delle righe locali
 * @param end_row Indice di fine delle righe locali
 */
void update_grid(double* U, double* U_new, int n, double h, double (*f)(double, double), int rank, int size, int start_row, int end_row) {

    // Communicate boundary rows with neighboring processes
    if (rank > 0) { // send to previous rank
        MPI_Send(U + n * start_row, n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(U + (start_row - 1) * n, n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank < size - 1) { // send to next rank
        MPI_Send(U + end_row * n, n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(U + (end_row + 1) * n, n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank > 0 && rank < size - 1) {
        //@note avoit pow for small integer exponents, use multiplication: pow(h,2) -> h*h
        // Even better, precompute h*h: it is a constant inside the loop. double h2 = h*h;
        // Also, if you use a local variable where you precompute the values of f at the node, you avoid multiple calls to the function
        #pragma omp parallel for
        for (int i = start_row; i < end_row; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                U_new[i * n + j] = (U[(i - 1) * n + j] + U[(i + 1) * n + j] + U[i * n + j - 1] + U[i * n + j + 1] + pow(h, 2) * f(static_cast<double>(i) / (n - 1), static_cast<double>(j) / (n - 1))) / 4.0;
            }
        }
    }
    // Caso per rank == 0
    else if (rank == 0) {
        #pragma omp parallel for
        for (int i = start_row + 1; i <= end_row; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                U_new[i * n + j] = (U[(i - 1) * n + j] + U[(i + 1) * n + j] + U[i * n + j - 1] + U[i * n + j + 1] + pow(h, 2) * f(static_cast<double>(i) / (n - 1), static_cast<double>(j) / (n - 1))) / 4.0;
            }
        }
    }
    // Caso per rank == size - 1
    else if (rank == size - 1) {
        
        #pragma omp parallel for
        for (int i = start_row; i <= end_row - 1; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                U_new[i * n + j] = (U[(i - 1) * n + j] + U[(i + 1) * n + j] + U[i * n + j - 1] + U[i * n + j + 1] + pow(h, 2) * f(static_cast<double>(i) / (n - 1), static_cast<double>(j) / (n - 1))) / 4.0;
            }
        }
    }
}
    
/**
 * @brief Controlla la convergenza della soluzione.
 * 
 * @param U Puntatore all'array che contiene i valori della griglia
 * @param U_new Puntatore all'array per i nuovi valori della griglia
 * @param n Dimensione della griglia
 * @param tol Tolleranza per la convergenza
 * @param rank Rank del processo MPI
 * @param size Numero totale di processi MPI
 * @param h Passo della griglia
 * @param start_row Indice di inizio delle righe locali
 * @param end_row Indice di fine delle righe locali
 * @return true se la soluzione è convergente, false altrimenti
 */
bool check_convergence(double* U, double* U_new, int n, double tol, int rank, int size, double h, int start_row, int end_row) {
    

    double local_sum = 0.0;

    #pragma omp parallel for reduction(+:local_sum)
    for (int i = start_row; i < end_row; ++i) {
        for (int j = 1; j < n - 1; ++j) {
            double diff = U_new[(i-start_row) * n + j] - U[(i-start_row) * n + j];
            local_sum += diff * diff;
        }
    }

    double global_sum;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    double error = sqrt(global_sum) * h;

    // std::cout << error << std::endl;
    return (error < tol);
}

/**
 * @brief Funzione di forzamento per l'equazione di Laplace.
 * 
 * @param x Coordinate x
 * @param y Coordinate y
 * @return Valore della funzione di forzamento
 */
double f_function(double x, double y) {
    return 8 * M_PI * M_PI * sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

/**
 * @brief Esporta i risultati della simulazione in formato VTK.
 * 
 * @param U Puntatore all'array che contiene i valori della griglia
 * @param n Dimensione della griglia
 * @param rank Rank del processo MPI
 * @param size Numero totale di processi MPI
 * @param num_threads Numero di thread OpenMP
 * @param from_script Se il codice è stato lanciato da uno script
 */
void export_to_vtk(double* U, int n, int rank, int size, int num_threads, bool from_script) {
    char filename[64];
    if (from_script){
    sprintf(filename, "test/data/solution_%d_%d_%d.vtk", n ,size, num_threads);
    } else {
        sprintf(filename, "solution_%d_%d_%d.vtk", n ,size, num_threads);
    }
    std::ofstream vtk_file;
    vtk_file.open(filename);

    vtk_file << "# vtk DataFile Version 2.0\n";
    vtk_file << "Laplace Solver Output\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET STRUCTURED_GRID\n";
    vtk_file << "DIMENSIONS " << n << " " << n << " 1\n";
    vtk_file << "POINTS " << n * n << " float\n";

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            vtk_file << i * 1.0 / (n - 1) << " " << j * 1.0 / (n - 1) << " 0\n";
        }
    }

    vtk_file << "POINT_DATA " << n * n << "\n";
    vtk_file << "SCALARS u float 1\n";
    vtk_file << "LOOKUP_TABLE default\n";

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            vtk_file << U[j * n + i] << "\n";
        }
    }

    vtk_file.close();
}



