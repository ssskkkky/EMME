#ifndef SOLVER_H
#define SOLVER_H
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

// #include <chrono>
#include <complex>
#include <iostream>
#include <vector>

#include "DedicatedThreadPool.h"
#include "Grid.h"
#include "Matrix.h"
#include "Parameters.h"
#include "lapack.h"

std::pair<std::complex<double>, Matrix<std::complex<double>>>
NewtonTraceIteration(std::complex<double> lambda, double tol);

std::pair<std::complex<double>, Matrix<std::complex<double>>>
NewtonTraceIterationSecantMethod(std::complex<double> lambda,
                                 const double& tol,
                                 Parameters& para,
                                 const Matrix<double>& coeff_matrix,
                                 const Grid<double>& grid_info,
                                 const int&);

template <typename T>

Matrix<T> NullSpace(const Matrix<T>& input, const double tol) {
    auto A = input;
    // Check if the matrix is square
    if (A.getRows() != A.getCols()) {
        throw std::invalid_argument("Input matrix must be square.");
    }

    // Perform Singular Value Decomposition (SVD)
    // You'll need an external library like Eigen or LAPACK for SVD
    // Replace this placeholder with your preferred SVD implementation
    // which returns the singular values (S), left singular vectors (U),
    // and right singular vectors (V)
    std::vector<double> S(A.getRows());  // why S is real??
    Matrix<T> U(A.getRows(), A.getRows());

    Matrix<T> VT(A.getCols(), A.getCols());

    const char* jobu = "None";
    const char* jobvt = "All";
    const lapack_int dimm = A.getRows();
    const lapack_int dimn = A.getCols();

    lapack_int work_length = dimm;
    // lapack_int optimal_work_length{};
    lapack_int lwork = 5 * dimm;
    std::vector<std::complex<double>> work(lwork);
    std::vector<double> rwork(5 * work_length);

    lapack_int info{};

    LAPACK_zgesvd(jobu, jobvt, &dimm, &dimn, A.data(), &dimm, S.data(),
                  U.data(), &dimm, VT.data(), &dimm, work.data(), &lwork,
                  rwork.data(), &info);

    auto V = VT;

    // ... perform SVD on A and store results in S, U, V

    // Identify null space basis vectors and construct the null space matrix
    int null_space_dim = 0;  // Count the number of null space basis vectors
    for (unsigned i = 0; i < A.getCols(); ++i) {
        // Check for singular values close to zero (tolerance approach)
        if (std::abs(S[i]) < tol) { null_space_dim++; }
    }

    // Create a result matrix to store the null space basis vectors as columns
    Matrix<T> null_space(A.getRows(), null_space_dim);
    int col_index = 0;
    for (unsigned int i = 0; i < A.getCols(); ++i) {
        // Check for singular values close to zero (tolerance approach)
        if (std::abs(S[i]) < tol) {
            // Extract the corresponding right singular vector and store in null
            // space matrix
            null_space.setCol(col_index, V.getCol(i));
            col_index++;
        }
    }
    return null_space;
};  // this is not tested

// Function representing the nonlinear eigenvalue problem (NLEP)
// F(lambda, x) = 0
template <typename Func, typename Funcb>
Matrix<std::complex<double>> F(const std::complex<double>& tau,
                               const std::complex<double>& beta_e,
                               const std::complex<double>& lambda,
                               const Func& kappa_f_tau_all,
                               const Funcb& bi,
                               const Matrix<double>& coeff_matrix,
                               const Grid<double>& grid_info) {
    // Implement the NLEP function here
    // This function should return a vector representing the residual (F(lambda,
    // x)) for a given eigenvalue (lambda) and eigenvector (x)
    // using namespace std::chrono;

    Matrix<std::complex<double>> quadrature_matrix(2 * grid_info.npoints,
                                                   2 * grid_info.npoints);
    Matrix<std::complex<double>> f_lambda(grid_info.npoints, grid_info.npoints);

    auto& thread_pool = DedicatedThreadPool<void>::get_instance();
    std::vector<std::future<void>> res;

    for (unsigned int i = 0; i < grid_info.npoints; i++) {
        for (unsigned int j = i; j < grid_info.npoints; j++) {
            if (i == j) {
                quadrature_matrix(i, j) = (1.0 + 1.0 / tau);
                quadrature_matrix(i, j + grid_info.npoints) = 0.0;
                quadrature_matrix(i + grid_info.npoints, j) = 0.0;
                quadrature_matrix(i + grid_info.npoints,
                                  j + grid_info.npoints) =
                    (2.0 * tau) / beta_e * bi(grid_info.grid[i]);
            } else {
                res.push_back(thread_pool.queue_task([&, i, j]() {
                    quadrature_matrix(i, j) =
                        -kappa_f_tau_all(0, grid_info.grid[i],
                                         grid_info.grid[j], lambda) *
                        coeff_matrix(i, j) * grid_info.dx;
                    quadrature_matrix(i, j + grid_info.npoints) =
                        kappa_f_tau_all(1, grid_info.grid[i], grid_info.grid[j],
                                        lambda) *
                        grid_info.dx;

                    // quadrature_matrix(i + grid_info.npoints, j) =
                    //     -kappa_f_tau_all(1, grid_info.grid[i],
                    //                      grid_info.grid[j], lambda) *
                    //     grid_info.dx;
                    quadrature_matrix(i + grid_info.npoints,
                                      j + grid_info.npoints) =
                        kappa_f_tau_all(2, grid_info.grid[i], grid_info.grid[j],
                                        lambda) *
                        grid_info.dx;

                    quadrature_matrix(j, i) = quadrature_matrix(i, j);

                    quadrature_matrix(j, i + grid_info.npoints) =
                        -quadrature_matrix(i, j + grid_info.npoints);
                    // quadrature_matrix(i + grid_info.npoints, j) =
                    //     -kappa_f_tau_all(1, grid_info.grid[i],
                    //                      grid_info.grid[j], lambda) *
                    //     grid_info.dx;
                    quadrature_matrix(j + grid_info.npoints,
                                      i + grid_info.npoints) =
                        quadrature_matrix(i + grid_info.npoints,
                                          j + grid_info.npoints);

                    quadrature_matrix(i + grid_info.npoints, j) =
                        quadrature_matrix(j, i + grid_info.npoints);

                    quadrature_matrix(j + grid_info.npoints, i) =
                        quadrature_matrix(i, j + grid_info.npoints);
                }));
            }
        }
    }
    for (auto& f : res) { f.get(); }

    return quadrature_matrix;
}

#endif  // SOLVER_H
