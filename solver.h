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
#include "aligned-allocator.h"
#include "lapack.h"

using value_type = std::complex<double>;
using matrix_type = Matrix<value_type, util::AlignedAllocator<value_type>>;

std::pair<value_type, matrix_type> NewtonTraceIteration(value_type lambda,
                                                        double tol);

std::pair<value_type, matrix_type> NewtonTraceIterationSecantMethod(
    value_type lambda,
    const double& tol,
    Parameters& para,
    const Matrix<double>& coeff_matrix,
    const Grid<double>& grid_info,
    const int&);

matrix_type NullSpace(const matrix_type& input, const double tol);

// Function representing the nonlinear eigenvalue problem (NLEP)
// F(lambda, x) = 0
template <typename FuncA, typename FuncB>
void F(const value_type& tau,
       const value_type& beta_e,
       const value_type& lambda,
       const FuncA& kappa_f_tau_all,
       const FuncB& bi,
       const Matrix<double>& coeff_matrix,
       const Grid<double>& grid_info,
       matrix_type& mat) {
    // Implement the NLEP function here
    // This function should return a vector representing the residual
    // (F(lambda, x)) for a given eigenvalue (lambda) and eigenvector (x)
    // using namespace std::chrono;
#ifdef EMME_DEBUG
    if (mat.getRows() != 2 * grid_info.npoints ||
        mat.getCols() != 2 * grid_info.npoints) {
        throw std::runtime_error("Matrix dimension and grid length mismatch.");
    }
#endif

    auto& thread_pool = DedicatedThreadPool<void>::get_instance();
    std::vector<std::future<void>> res;

    for (unsigned int i = 0; i < grid_info.npoints; i++) {
        for (unsigned int j = i; j < grid_info.npoints; j++) {
            if (i == j) {
                mat(i, j) = (1.0 + 1.0 / tau);
                mat(i, j + grid_info.npoints) = 0.0;
                mat(i + grid_info.npoints, j) = 0.0;
                mat(i + grid_info.npoints, j + grid_info.npoints) =
                    (2.0 * tau) / beta_e * bi(grid_info.grid[i]);
            } else {
                res.push_back(thread_pool.queue_task([&, i, j]() {
                    mat(i, j) = -kappa_f_tau_all(0, grid_info.grid[i],
                                                 grid_info.grid[j], lambda) *
                                coeff_matrix(i, j) * grid_info.dx;
                    mat(i, j + grid_info.npoints) =
                        kappa_f_tau_all(1, grid_info.grid[i], grid_info.grid[j],
                                        lambda) *
                        grid_info.dx;

                    mat(i + grid_info.npoints, j + grid_info.npoints) =
                        kappa_f_tau_all(2, grid_info.grid[i], grid_info.grid[j],
                                        lambda) *
                        grid_info.dx;

                    mat(j, i) = mat(i, j);

                    mat(j, i + grid_info.npoints) =
                        -mat(i, j + grid_info.npoints);
                    mat(j + grid_info.npoints, i + grid_info.npoints) =
                        mat(i + grid_info.npoints, j + grid_info.npoints);

                    mat(i + grid_info.npoints, j) =
                        mat(j, i + grid_info.npoints);

                    mat(j + grid_info.npoints, i) =
                        mat(i, j + grid_info.npoints);
                }));
            }
        }
    }
    for (auto& f : res) { f.get(); }
}

template <typename FuncA, typename FuncB>
matrix_type F(const value_type& tau,
              const value_type& beta_e,
              const value_type& lambda,
              const FuncA& kappa_f_tau_all,
              const FuncB& bi,
              const Matrix<double>& coeff_matrix,
              const Grid<double>& grid_info) {
    matrix_type quadrature_matrix(2 * grid_info.npoints, 2 * grid_info.npoints);
    F(tau, beta_e, lambda, kappa_f_tau_all, bi, coeff_matrix, grid_info,
      quadrature_matrix);
    return quadrature_matrix;
}

#endif  // SOLVER_H
