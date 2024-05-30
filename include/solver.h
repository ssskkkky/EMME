#ifndef SOLVER_H
#define SOLVER_H
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

// #include <chrono>
#include <complex>
#include <vector>

#include "DedicatedThreadPool.h"
#include "Grid.h"
#include "Matrix.h"
#include "Parameters.h"
#include "Timer.h"
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

template <typename T>
class EigenSolver {
   public:
    using value_type = typename T::value_type;
    using matrix_type = T;
    // EigenSolver(const Parameters& para_input,
    //             value_type eigen_init,
    //             const Matrix<double>& coeff_matrix_input,
    //             const Grid<double>& grid_info_input){};
    // matrix_type matrixAssembler(Parameters, value_type, Matrix<double>){};
    void matrixDerivativeSecantAssembler() {
        eigen_matrix_derivative =
            (eigen_matrix - eigen_matrix_old) / d_eigen_value;
    };
    std::vector<value_type> nullSpace() {
        auto A = eigen_matrix;
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
        matrix_type U(A.getRows(), A.getRows());

        matrix_type VT(A.getCols(), A.getCols());

        const char* jobz = "All";
        const lapack_int dimm = A.getRows();
        const lapack_int dimn = A.getCols();

        lapack_int work_length = dimm;
        // lapack_int optimal_work_length{};
        lapack_int lwork = 2 * work_length * work_length;
        std::vector<value_type> work(lwork);
        std::vector<double> rwork(50 * work_length * work_length);
        std::vector<int> iwork(8 * work_length);

        lapack_int info{};

        // const char* jobu = "None";
        // const char* jobvt = "All";
        // LAPACK_zgesvd(jobu, jobvt, &dimm, &dimn, A.data(), &dimm, S.data(),
        //               U.data(), &dimm, VT.data(), &dimm, work.data(), &lwork,
        //               rwork.data(), &info);

        // zgesdd is much faster than zgesvd
        LAPACK_zgesdd(jobz, &dimm, &dimn, A.data(), &dimm, S.data(), U.data(),
                      &dimm, VT.data(), &dimm, work.data(), &lwork,
                      rwork.data(), iwork.data(), &info);

        auto V = VT;

        // ... perform SVD on A and store results in S, U, V

        // Identify null space basis vectors and construct the null space matrix
        int null_space_dim = 0;  // Count the number of null space basis vectors
        for (unsigned i = 0; i < A.getCols(); ++i) {
            // Check for singular values close to zero (tolerance approach)
            if (std::abs(S[i]) < null_space_tol) { null_space_dim++; }
        }

        // Create a result matrix to store the null space basis vectors as
        // columns
        matrix_type null_space(A.getRows(), null_space_dim);
        int col_index = 0;
        for (unsigned int i = 0; i < A.getCols(); ++i) {
            // Check for singular values close to zero (tolerance approach)
            if (std::abs(S[i]) < null_space_tol) {
                // Extract the corresponding right singular vector and store in
                // null space matrix
                null_space.setCol(col_index, V.getCol(i));
                col_index++;
            }
        }
        return V.getCol(V.getCols() - 1);
    };
    void newtonTraceSecantIteration() {
        eigen_matrix_old = eigen_matrix;
        const char* upper = "Upper";
        const lapack_int dim = eigen_matrix.getRows();
        lapack_int work_length = dim * dim;
        lapack_int optimal_work_length{};
        std::vector<value_type> work(work_length);
        std::vector<lapack_int> ipiv(dim);
        lapack_int info{};

        if (optimal_work_length > work_length) {
            work_length = optimal_work_length;
            work.resize(work_length);
        }

        Timer::get_timer().start_timing("linear solver");
        LAPACK_zsysv(upper, &dim, &dim, eigen_matrix.data(), &dim, ipiv.data(),
                     eigen_matrix_derivative.data(), &dim, work.data(),
                     &work_length, &info);
        Timer::get_timer().pause_timing("linear solver");
        d_eigen_value = -1.0 / eigen_matrix_derivative.trace();
        eigen_value += d_eigen_value;

        // if (info != 0) {
        //     if (info < 0) {
        //         std::cout << "the " << -info
        //                   << "-th argument had an illegal value";
        //     } else {
        //         std::cout << "The factorization has been completed, but the "
        //                      "block diagonal matrix D is exactly singular at
        //                      "
        //                   << i << ", so the solution could not be computed.";
        //     }
        //     throw std::runtime_error("Linear solve failed");
        // };
        optimal_work_length = work[0].real();
        Timer::get_timer().start_timing("integration");
        matrixAssembler(eigen_matrix);
        Timer::get_timer().pause_timing("integration");
        matrixDerivativeSecantAssembler();
    }
    const Parameters& para;
    value_type eigen_value;
    value_type d_eigen_value;
    double null_space_tol;
    const Matrix<double>& coeff_matrix;
    const Grid<double>& grid_info;
    matrix_type eigen_matrix;
    matrix_type eigen_matrix_old;
    matrix_type eigen_matrix_derivative;

    EigenSolver(const Parameters& para_input,
                value_type eigen_init,
                const Matrix<double>& coeff_matrix_input,
                const Grid<double>& grid_info_input)
        : para(para_input),
          eigen_value(0.99 * eigen_init),
          d_eigen_value(0.01 * eigen_init),
          null_space_tol(1e-1),
          coeff_matrix(coeff_matrix_input),
          grid_info(grid_info_input),
          eigen_matrix(grid_info.npoints * 2, 2 * grid_info.npoints),
          eigen_matrix_old(grid_info.npoints * 2, 2 * grid_info.npoints),
          eigen_matrix_derivative(grid_info.npoints * 2,
                                  2 * grid_info.npoints) {
        matrixAssembler(eigen_matrix_old);
        eigen_value += d_eigen_value;
        matrixAssembler(eigen_matrix);
        matrixDerivativeSecantAssembler();
    }

    void matrixAssembler(matrix_type& mat) {
#ifdef EMME_DEBUG
        if (mat.getRows() != 2 * grid_info.npoints ||
            mat.getCols() != 2 * grid_info.npoints) {
            throw std::runtime_error(
                "Matrix dimension and grid length mismatch.");
        }
#endif

        static auto kappa_f_tau_all = [&](unsigned i, double eta, double eta_p,
                                          std::complex<double> omega) {
            return para.kappa_f_tau(i, eta, eta_p, omega) +
                   para.kappa_f_tau_e(i, eta, eta_p, omega);
        };

#ifdef MULTI_THREAD
        auto& thread_pool = DedicatedThreadPool<void>::get_instance();

        std::vector<std::future<void>> res;
#endif

        for (unsigned int i = 0; i < grid_info.npoints; i++) {
            for (unsigned int j = i; j < grid_info.npoints; j++) {
                if (i == j) {
                    mat(i, j) = (1.0 + 1.0 / para.tau);
                    mat(i, j + grid_info.npoints) = 0.0;
                    mat(i + grid_info.npoints, j) = 0.0;
                    mat(i + grid_info.npoints, j + grid_info.npoints) =
                        (2.0 * para.tau) / para.beta_e *
                        para.bi(grid_info.grid[i]);
                } else {
#ifdef MULTI_THREAD
                    res.push_back(thread_pool.queue_task([&, i, j]() {
#endif
                        mat(i, j) =
                            -kappa_f_tau_all(0, grid_info.grid[i],
                                             grid_info.grid[j], eigen_value) *
                            coeff_matrix(i, j) * grid_info.dx;
                        mat(i, j + grid_info.npoints) =
                            kappa_f_tau_all(1, grid_info.grid[i],
                                            grid_info.grid[j], eigen_value) *
                            grid_info.dx;

                        mat(i + grid_info.npoints, j + grid_info.npoints) =
                            kappa_f_tau_all(2, grid_info.grid[i],
                                            grid_info.grid[j], eigen_value) *
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
#ifdef MULTI_THREAD
                    }));
#endif
                }
            }
        }
#ifdef MULTI_THREAD
        for (auto& f : res) { f.get(); }
#endif
    }
};
#endif  // SOLVER_H
