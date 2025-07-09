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
#include <fstream>
#include <nlohmann/json.hpp>

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

    // Utility function to export a matrix to a JSON file
    // void export_matrix_to_json(const matrix_type& mat, const std::string& label, const std::string& step, const std::string& filename = "qr_debug.json") {
    //     nlohmann::json j;
    //     j["step"] = step;
    //     j["label"] = label;
    //     j["rows"] = mat.getRows();
    //     j["cols"] = mat.getCols();
    //     j["data"] = nlohmann::json::array();
    //     for (unsigned i = 0; i < mat.getRows(); ++i) {
    //         nlohmann::json row = nlohmann::json::array();
    //         for (unsigned jcol = 0; jcol < mat.getCols(); ++jcol) {
    //             row.push_back({mat(i, jcol).real(), mat(i, jcol).imag()});
    //         }
    //         j["data"].push_back(row);
    //     }
    //     std::ofstream ofs(filename, std::ios::app);
    //     ofs << j.dump() << std::endl;
    // }

    // Utility function to export a vector to a JSON file
    // void export_vector_to_json(const std::vector<value_type>& vec, const std::string& label, const std::string& step, const std::string& filename = "qr_debug.json") {
    //     nlohmann::json j;
    //     j["step"] = step;
    //     j["label"] = label;
    //     j["size"] = vec.size();
    //     j["data"] = nlohmann::json::array();
    //     for (const auto& v : vec) {
    //         j["data"].push_back({v.real(), v.imag()});
    //     }
    //     std::ofstream ofs(filename, std::ios::app);
    //     ofs << j.dump() << std::endl;
    // }

    // inline void export_scalar_to_json(const value_type& val, const std::string& label, const std::string& step, const std::string& filename = "qr_debug.json") {
    //     nlohmann::json j;
    //     j["label"] = label;
    //     j["step"] = step;
    //     j["value"] = {val.real(), val.imag()};
    //     std::ofstream file(filename, std::ios::app);
    //     file << j.dump() << std::endl;
    // }

    void newtonQRSecantIteration() {
        eigen_matrix_old = eigen_matrix;
        const lapack_int dim = eigen_matrix.getRows();
        lapack_int work_length = dim * dim;
        lapack_int optimal_work_length{};
        std::vector<value_type> work(work_length);
        std::vector<lapack_int> ipiv(dim);
        lapack_int info{};
        std::vector<value_type> tau(dim);         // Added declaration for tau
        lapack_int lwork = work_length;           // Added declaration for lwork
        std::vector<value_type> R_last_col(dim);  // Add this before first use
        std::vector<lapack_int> jpvt(
            dim);  // Permutation vector for column pivoting
        std::vector<double> rwork(2 * dim);  // Real work array for zgeqp3

        // Initialize permutation vector (no initial pivoting)
        for (int i = 0; i < dim; ++i) { jpvt[i] = 0; }

        if (optimal_work_length > work_length) {
            work_length = optimal_work_length;
            work.resize(work_length);
        }

        Timer::get_timer().start_timing("linear solver");
        // All export_matrix_to_json, export_vector_to_json, and
        // export_scalar_to_json calls are removed from this function and
        // related debug output functions. No debug files will be created or
        // written.
#ifdef EMME_MKL
        matrix_type eigen_matrix_col_major(dim, dim);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                eigen_matrix_col_major(j, i) = eigen_matrix(i, j);

        zgeqp3(&dim, &dim, eigen_matrix_col_major.data(), &dim, jpvt.data(),
               tau.data(), work.data(), &lwork, rwork.data(), &info);
#else
        matrix_type eigen_matrix_col_major(dim, dim);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                eigen_matrix_col_major(j, i) = eigen_matrix(i, j);

        LAPACK_zgeqp3(&dim, &dim, eigen_matrix_col_major.data(), &dim,
                      jpvt.data(), tau.data(), work.data(), &lwork,
                      rwork.data(), &info);
#endif
        if (info != 0)
            throw std::runtime_error("QR factorization with pivoting failed");

        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                eigen_matrix(i, j) = eigen_matrix_col_major(j, i);

        // After QR factorization, print diagonal of eigen_matrix
        // std::cout << "[DEBUG] Diagonal of eigen_matrix after QR
        // factorization:" << std::endl; for (int i = 0; i < dim; ++i) {
        //     std::cout << i << ": " << eigen_matrix(i, i) << std::endl;
        // }

        // Print permutation vector for debugging
        // std::cout << "[DEBUG] Column permutation vector jpvt:" << std::endl;
        // for (int i = 0; i < dim; ++i) {
        //     std::cout << i << ": " << jpvt[i] << std::endl;
        // }

        // Find which column corresponds to the last column in the original
        // matrix
        int last_col_permuted = -1;
        for (int i = 0; i < dim; ++i) {
            if (jpvt[i] - 1 == dim - 1) {  // jpvt is 1-indexed
                last_col_permuted = i;
                break;
            }
        }
        if (last_col_permuted == -1) {
            throw std::runtime_error("Could not find permuted last column");
        }

        matrix_type R_mat_except_last_col(dim - 1, dim - 1);
        for (int col = 0; col < dim - 1; ++col) {
            for (int row = 0; row <= col; ++row) {
                R_mat_except_last_col(row, col) = eigen_matrix(row, col);
            }
        }
        for (int row = 0; row <= dim - 2; ++row) {
            R_last_col[row] = eigen_matrix(row, dim - 1);
        }
        // export_matrix_to_json(R_mat_except_last_col, "R_mat_except_last_col",
        // "before_triangular_solve");
        matrix_type R_last_col_mat_before(dim, 1);
        for (int i = 0; i < dim - 1; ++i) {
            R_last_col_mat_before(i, 0) = R_last_col[i];
        }
        R_last_col_mat_before(dim - 1, 0) = -1.0;
        // export_matrix_to_json(R_last_col_mat_before, "R_last_col",
        // "before_triangular_solve");

        const char uplo = 'U';
        const char trans = 'N';
        const char diag = 'N';
        lapack_int one = 1;
        lapack_int n_rhs = 1;
        lapack_int r_dim = dim - 1;
        matrix_type R_mat_except_last_col_col_major(r_dim, r_dim);
        for (int i = 0; i < r_dim; ++i)
            for (int j = 0; j < r_dim; ++j)
                R_mat_except_last_col_col_major(j, i) =
                    R_mat_except_last_col(i, j);
        LAPACK_ztrtrs(&uplo, &trans, &diag, &r_dim, &n_rhs,
                      R_mat_except_last_col_col_major.data(), &r_dim,
                      R_last_col.data(), &r_dim, &info);
        if (info != 0) {
            if (info < 0) {
                throw std::runtime_error("第 " + std::to_string(-info) +
                                         " 个参数非法");
            } else {
                throw std::runtime_error("矩阵第 " + std::to_string(info) +
                                         " 个对角线元素为零，无法求解");
            }
        }
        // export_matrix_to_json(R_mat_except_last_col, "R_mat_except_last_col",
        // "after_triangular_solve");
        matrix_type R_last_col_mat_after(dim, 1);
        for (int i = 0; i < dim; ++i) {
            R_last_col_mat_after(i, 0) = R_last_col[i];
        }
        // export_matrix_to_json(R_last_col_mat_after, "R_last_col_mat_after",
        // "after_triangular_solve");

        const char side = 'L';  // 应用 Q 到左侧
        // Remove duplicate declaration of trans after LAPACK_ztrtrs

        matrix_type q_last(dim, 1);  // Use correct constructor
        for (int i = 0; i < dim; ++i) q_last(i, 0) = 0.0;  // Zero-initialize
        q_last(last_col_permuted, 0) =
            1.0;  // Set element corresponding to permuted last column

        // 调用 zunmqr 应用 Q^H (共轭转置)
        LAPACK_zunmqr(&side, &trans, &dim, &one, &dim, eigen_matrix.data(),
                      &dim, tau.data(), q_last.data(), &dim, work.data(),
                      &lwork, &info);
        if (info != 0) throw std::runtime_error("Q application failed");

        Timer::get_timer().pause_timing("linear solver");
        // --- Correction: Extract R and Q for update ---
        // R is upper triangular in eigen_matrix (row-major after transpose
        // back) Q is not explicitly formed, but we can use LAPACK's orgqr/unmqr
        // to get Q We'll use the last row of Q (after permutation) and the last
        // element of R

        // 1. Extract R (upper triangle of eigen_matrix)
        matrix_type R_mat(dim, dim);
        for (int i = 0; i < dim; ++i)
            for (int j = i; j < dim; ++j) R_mat(i, j) = eigen_matrix(i, j);

        // export_matrix_to_json(R_mat, "R_mat", "after_R_formation");

        // 2. Form Q explicitly (column-major for LAPACK)
        matrix_type Q_mat(dim, dim);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j) Q_mat(j, i) = eigen_matrix(j, i);
        // Use LAPACK_zungqr to form Q from the output of zgeqp3
        lapack_int info_q = 0;
#ifdef EMME_MKL
        matrix_type Q_mat_col_major(dim, dim);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                Q_mat_col_major(j, i) = eigen_matrix(i, j);

        zungqr(&dim, &dim, &dim, Q_mat_col_major.data(), &dim, tau.data(),
               work.data(), &lwork, &info_q);
#else
        matrix_type Q_mat_col_major(dim, dim);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                Q_mat_col_major(j, i) = eigen_matrix(i, j);

        LAPACK_zungqr(&dim, &dim, &dim, Q_mat_col_major.data(), &dim,
                      tau.data(), work.data(), &lwork, &info_q);
#endif

        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j) Q_mat(i, j) = Q_mat_col_major(j, i);

        if (info_q != 0) throw std::runtime_error("Q formation failed");
        // Export Q_mat as its conjugate transpose for compatibility with
        // Mathematica
        matrix_type Q_mat_conjT(dim, dim);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                Q_mat_conjT(j, i) = std::conj(Q_mat(j, i));
        // export_matrix_to_json(Q_mat_conjT, "Q_mat", "after_Q_formation");

        Q_mat = Q_mat_conjT;
        // export_matrix_to_json(Q_mat, "Q_mat", "after_Q_formation");

        // 3. Numerator: last element of permuted R
        value_type R_last = R_mat(dim - 1, dim - 1);
        // export_scalar_to_json(R_last, "R_last", "after_R_last_computation");
        // 4. Construct the vector as in Mathematica: v_full = Append[-p, 1.0]
        std::vector<value_type> v_full(dim);
        for (int i = 0; i < dim - 1; ++i) { v_full[i] = -R_last_col[i]; }
        v_full[dim - 1] = 1.0;
        // export_vector_to_json(v_full, "v_full", "after_v_full_construction");
        // 5. Permute v_full according to jpvt
        std::vector<value_type> v_full_permuted(dim);
        for (int i = 0; i < dim; ++i) {
            v_full_permuted[jpvt[i] - 1] = v_full[i];  // jpvt is 1-based
        }
        // export_vector_to_json(v_full_permuted, "v_full_permuted",
        // "after_v_full_permutation");
        // 6. Denominator: last row of permuted Q dotted with
        // (eigen_matrix_derivative * v_full_permuted)
        // export_matrix_to_json(eigen_matrix_derivative,
        // "eigen_matrix_derivative", "before_derivative_update");
        std::vector<value_type> temp_vec(dim, 0.0);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                temp_vec[i] +=
                    eigen_matrix_derivative(i, j) * v_full_permuted[j];
            }
        }
        // export_vector_to_json(temp_vec, "temp_vec",
        // "after_temp_vec_computation");
        // export_matrix_to_json(eigen_matrix_derivative,
        // "eigen_matrix_derivative", "after_derivative_update");

        // Calculate Q_dot as the dot product of the last column of Q with
        // temp_vec
        value_type Q_dot = 0.0;
        for (int i = 0; i < dim; ++i) {
            Q_dot += Q_mat(i, dim - 1) * temp_vec[i];
        }

        // export_scalar_to_json(R_last, "R_last", "after_R_last_computation");
        // export_scalar_to_json(Q_dot, "Q_dot", "after_Q_dot_computation");

        d_eigen_value = -R_last / Q_dot;
        // export_scalar_to_json(d_eigen_value, "d_eigen_value",
        // "after_d_eigen_value_computation");
        eigen_value += d_eigen_value;
        // export_matrix_to_json(eigen_matrix, "eigen_matrix",
        // "after_eigenvalue_update");

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
