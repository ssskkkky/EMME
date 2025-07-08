#ifndef SOLVER_H
#define SOLVER_H
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#ifdef EMME_MKL
#define MKL_Complex16 lapack_complex_double
#define lapack_int MKL_INT  // This do work
#endif

// #include <chrono>
#include <complex>
#include <iostream>
#include <vector>

#include "DedicatedThreadPool.h"
#include "Grid.h"
#include "Matrix.h"
#include "Parameters.h"
#include "Timer.h"
#include "aligned-allocator.h"
#ifdef EMME_MKL
#include "mkl_lapack.h"
#else
#include "lapack.h"
#endif

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
#ifdef EMME_MKL
        zgesdd(jobz, &dimm, &dimn, A.data(), &dimm, S.data(), U.data(), &dimm,
               VT.data(), &dimm, work.data(), &lwork, rwork.data(),
               iwork.data(), &info);

#else
        LAPACK_zgesdd(jobz, &dimm, &dimn, A.data(), &dimm, S.data(), U.data(),
                      &dimm, VT.data(), &dimm, work.data(), &lwork,
                      rwork.data(), iwork.data(), &info);
#endif
        // lapack is coloum major order but we are Row major order, VT is not V
        auto kernel_vector = VT.getCol(VT.getCols() - 1);

        for (auto& ele : kernel_vector) { ele = std::conj(ele); }

        return kernel_vector;
        ;
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
#ifdef EMME_MKL
        zsysv(upper, &dim, &dim, eigen_matrix.data(), &dim, ipiv.data(),
              eigen_matrix_derivative.data(), &dim, work.data(), &work_length,
              &info);
#else
        LAPACK_zsysv(upper, &dim, &dim, eigen_matrix.data(), &dim, ipiv.data(),
                     eigen_matrix_derivative.data(), &dim, work.data(),
                     &work_length, &info);
#endif
        Timer::get_timer().pause_timing("linear solver");
        d_eigen_value = -1.0 / eigen_matrix_derivative.trace();
        eigen_value += d_eigen_value;

        if (info != 0) {
            std::ostringstream oss;
            oss << "Linear solve failed. ";
            if (info < 0) {
                oss << "the " << -info << "-th argument had an illegal value";
            } else {
                oss << "The factorization has been completed, but the "
                    << "block diagonal matrix D is exactly singular at " << info
                    << ", so the solution could not be computed.";
            }
            throw std::runtime_error(oss.str());
        };

        optimal_work_length = work[0].real();
        Timer::get_timer().start_timing("integration");
        matrixAssembler(eigen_matrix);
        Timer::get_timer().pause_timing("integration");
        matrixDerivativeSecantAssembler();
    }

    void newtonQRSecantIteration() {
        eigen_matrix_old = eigen_matrix;
        const char* upper = "Upper";
        const lapack_int dim = eigen_matrix.getRows();
        lapack_int work_length = dim * dim;
        lapack_int optimal_work_length{};
        std::vector<value_type> work(work_length);
        std::vector<lapack_int> ipiv(dim);
        lapack_int info{};
        std::vector<value_type> tau(dim); // Added declaration for tau
        lapack_int lwork = work_length;  // Added declaration for lwork
        std::vector<value_type> R_last_col(dim); // Add this before first use

        if (optimal_work_length > work_length) {
            work_length = optimal_work_length;
            work.resize(work_length);
        }

        std::cout << "[DEBUG] Entered newtonQRSecantIteration" << std::endl;
        std::cout << "[DEBUG] About to print diagonal elements" << std::endl;
        std::cout << "Diagonal elements of eigen_matrix before QR factorization:\n";
        for (int i = 0; i < dim; ++i) {
            std::cout << i << ": " << eigen_matrix(i, i) << std::endl;
        }
        Timer::get_timer().start_timing("linear solver");
#ifdef EMME_MKL
        zgeqrf(&dim, &dim, eigen_matrix.data(), &dim, tau.data(), work.data(),
               &lwork, &info);
#else

        LAPACK_zgeqrf(&dim, &dim, eigen_matrix.data(), &dim, tau.data(),
                      work.data(), &lwork, &info);
#endif
        if (info != 0) throw std::runtime_error("QR factorization failed");

        // After QR factorization, print diagonal of eigen_matrix
        std::cout << "[DEBUG] Diagonal of eigen_matrix after QR factorization:" << std::endl;
        for (int i = 0; i < dim; ++i) {
            std::cout << i << ": " << eigen_matrix(i, i) << std::endl;
        }

        matrix_type R_mat_except_last_col(dim - 1, dim - 1);
        int n = dim, m = dim;
        for (int col = 0; col < dim - 2; ++col) {
            for (int row = 0; row <= col; ++row) {  // 上三角部分
                R_mat_except_last_col(row, col) =
                    eigen_matrix(row, col);  // Use operator() for matrix access
            }
        }

        for (int row = 0; row <= dim - 2; ++row) {  // 上三角部分
            R_last_col[row] =
                eigen_matrix(row, dim - 1);  // Use operator() for matrix access
        }
        // Do not set R_last_col[dim - 1] = -1.0; Keep its value as is

        // Before triangular solve, print diagonal of R_mat_except_last_col
        std::cout << "[DEBUG] Diagonal of R_mat_except_last_col before triangular solve:" << std::endl;
        for (int i = 0; i < dim - 1; ++i) {
            std::cout << i << ": " << R_mat_except_last_col(i, i) << std::endl;
        }

        const char uplo = 'U';   // 上三角矩阵
        const char trans = 'N';  // 不进行转置 (A*X = B)
        const char diag = 'N';   // 是否单位三角矩阵

        // 调用 LAPACK 函数 ztrtrs
        lapack_int one = 1;
        lapack_int n_rhs = 1;
        lapack_int r_dim = dim - 1;
        LAPACK_ztrtrs(&uplo, &trans, &diag, &r_dim, &n_rhs,
                      R_mat_except_last_col.data(), &r_dim, R_last_col.data(), &r_dim,
                      &info);

        // 错误处理
        if (info != 0) {
            std::cout << "[DEBUG] Error: info = " << info << ", R_mat_except_last_col(" << (info-1) << ", " << (info-1) << ") = " << R_mat_except_last_col(info-1, info-1) << std::endl;
            if (info < 0) {
                throw std::runtime_error("第 " + std::to_string(-info) +
                                         " 个参数非法");
            } else {
                throw std::runtime_error("矩阵第 " + std::to_string(info) +
                                         " 个对角线元素为零，无法求解");
            }
        }

        const char side = 'L';   // 应用 Q 到左侧
        // Remove duplicate declaration of trans after LAPACK_ztrtrs

        matrix_type q_last(dim, 1); // Use correct constructor
        for (int i = 0; i < dim; ++i) q_last(i, 0) = 0.0; // Zero-initialize
        q_last(dim-1, 0) = 1.0;  // Set last element

        // 调用 zunmqr 应用 Q^H (共轭转置)
        LAPACK_zunmqr(&side, &trans, &dim, &one, &dim, eigen_matrix.data(), &dim,
                      tau.data(), q_last.data(), &dim, work.data(), &lwork,
                      &info);
        if (info != 0) throw std::runtime_error("Q application failed");

        Timer::get_timer().pause_timing("linear solver");
        d_eigen_value = -1.0 / eigen_matrix_derivative.trace();

        // Remove duplicate declaration of trans
        // Fix d_eigen_value calculation
        // d_eigen_value = -eigen_matrix[dim - 1, dim - 1] /
        //                (eigen_matrix_derivative.(-1.0 * R_last_col));
        // Assuming eigen_matrix_derivative is a matrix and R_last_col is a vector,
        // and operator* is overloaded for matrix-vector multiplication:
        std::vector<value_type> temp_vec(dim, 0.0);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                temp_vec[i] += eigen_matrix_derivative(i, j) * R_last_col[j];
            }
        }
        d_eigen_value = -eigen_matrix(dim - 1, dim - 1) / temp_vec[dim - 1];
        eigen_value += d_eigen_value;

        if (info != 0) {
            std::ostringstream oss;
            oss << "Linear solve failed. ";
            if (info < 0) {
                oss << "the " << -info << "-th argument had an illegal value";
            } else {
                oss << "The factorization has been completed, but the "
                    << "block diagonal matrix D is exactly singular at " << info
                    << ", so the solution could not be computed.";
            }
            throw std::runtime_error(oss.str());
        };

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
    unsigned int dim;
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
          dim(std::fpclassify(para.beta_e) == FP_ZERO ? grid_info.npoints
                                                      : 2 * grid_info.npoints),
          eigen_matrix(dim, dim),
          eigen_matrix_old(dim, dim),
          eigen_matrix_derivative(dim, dim) {
        matrixAssembler(eigen_matrix_old);
        eigen_value += d_eigen_value;
        matrixAssembler(eigen_matrix);
        matrixDerivativeSecantAssembler();
    }

    void matrixAssembler(matrix_type& mat) {
#ifdef EMME_DEBUG

        if (mat.getRows() != dim || mat.getCols() != dim) {
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

        if (std::fpclassify(para.beta_e) == FP_ZERO) {
            for (unsigned int i = 0; i < grid_info.npoints; i++) {
                for (unsigned int j = i; j < grid_info.npoints; j++) {
                    if (i == j) {
                        mat(i, j) = (1.0 + 1.0 / para.tau);
                    } else {
#ifdef MULTI_THREAD
                        res.push_back(thread_pool.queue_task([&, i, j]() {
#endif
                            mat(i, j) = -kappa_f_tau_all(0, grid_info.grid[i],
                                                         grid_info.grid[j],
                                                         eigen_value) *
                                        coeff_matrix(i, j) * grid_info.dx;

                            mat(j, i) = mat(i, j);
#ifdef MULTI_THREAD
                        }));
#endif
                    }
                }
            }

        } else {
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
                            mat(i, j) = -kappa_f_tau_all(0, grid_info.grid[i],
                                                         grid_info.grid[j],
                                                         eigen_value) *
                                        coeff_matrix(i, j) * grid_info.dx;
                            mat(i, j + grid_info.npoints) =
                                kappa_f_tau_all(1, grid_info.grid[i],
                                                grid_info.grid[j],
                                                eigen_value) *
                                grid_info.dx;

                            mat(i + grid_info.npoints, j + grid_info.npoints) =
                                kappa_f_tau_all(2, grid_info.grid[i],
                                                grid_info.grid[j],
                                                eigen_value) *
                                grid_info.dx;

                            mat(j, i) = mat(i, j);

                            mat(j, i + grid_info.npoints) =
                                -mat(i, j + grid_info.npoints);
                            mat(j + grid_info.npoints, i + grid_info.npoints) =
                                mat(i + grid_info.npoints,
                                    j + grid_info.npoints);

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
        }
#ifdef MULTI_THREAD
        for (auto& f : res) { f.get(); }
#endif
    }
};
#endif  // SOLVER_H
