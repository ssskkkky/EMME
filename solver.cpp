#include "solver.h"
#include <complex>
#include <iostream>
#include <vector>
#include "Grid.h"
#include "Matrix.h"
#include "Parameters.h"
#include "functions.h"

// Jacobian of F(lambda)
Matrix<std::complex<double>> Jacobian(std::complex<double> lambda) {
    // Implement the Jacobian of F(lambda, x) here
    // This function should return a 2D vector representing the Jacobian matrix
    Matrix<std::complex<double>> j(2, 2);
    // ... (fill the Jacobian matrix with partial derivatives)
    j(0, 0) = 2.0 * std::complex<double>(0, 1) * lambda *
              std::exp(std::complex<double>(0, 1) * lambda * lambda);
    j(0, 1) = 0.0;
    j(1, 0) = 0.0;
    j(1, 1) = 0.0;

    return j;
}

// Solve the linear system J * dx = -F using LU decomposition

std::vector<std::complex<double>> SolveLinearSystem(
    const Matrix<std::complex<double>>& J,
    const std::vector<std::complex<double>>& b)

{
    auto Jv = J;
    auto bv = b;
    int n = J.getRows();

    // Allocate memory for L, U, perm (assuming they have the same size as J)
    Matrix<std::complex<double>> L(n, n);
    Matrix<std::complex<double>> U(n, n);
    std::vector<int> perm(n);

    for (int i = 0; i < n; ++i) { perm[i] = i; }

    // Perform LU decomposition (manual implementation)
    for (int i = 0; i < n; ++i) {
        // Partial pivoting (find the largest element in the current column
        // below the diagonal)
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(Jv(j, i)) > std::abs(Jv(pivot, i))) { pivot = j; }
        }

        // Swap rows if necessary (update permutation vector)
        if (pivot != i) {
            std::swap(perm[i], perm[pivot]);
            for (int j = 0; j < n; ++j) { std::swap(Jv(i, j), Jv(pivot, j)); }
            std::swap(bv[i], bv[pivot]);
        }

        // LU decomposition for the current column
        for (int j = i + 1; j < n; ++j) {
            std::complex<double> factor = Jv(j, i) / Jv(i, i);
            Jv(j, i) = factor;
            for (int k = i + 1; k < n; ++k) { Jv(j, k) -= factor * Jv(i, k); }
        }
    }

    // Extract L and U from the modified Jv matrix
    for (int i = 0; i < n; ++i) {
        L(i, i) = 1;
        for (int j = 0; j < i; ++j) { L(i, j) = Jv(i, j); }
        for (int j = i; j < n; ++j) { U(i, j) = Jv(i, j); }
    }

    // for (const auto& element : L) { std::cout << element << " "; }
    // std::cout << std::endl;
    // for (const auto& element : U) { std::cout << element << " "; }
    // std::cout << std::endl;

    // Forward substitution (solve Ly = b)
    std::vector<std::complex<double>> y(n);
    for (int i = 0; i < n; ++i) {
        std::complex<double> sum = 0.0;
        for (int j = 0; j < i; ++j) { sum += L(i, j) * y[j]; }
        //       for (int j = 0; j < i; ++j) { sum += L(i, j) * y[j]; }
        y[i] = bv[i] - sum;
    }

    // Backward substitution (solve Ux = y)
    std::vector<std::complex<double>> x(n);
    for (int i = n - 1; i >= 0; --i) {
        std::complex<double> sum = 0.0;
        for (int j = i + 1; j < n; ++j) { sum += U(i, j) * x[j]; }
        x[i] = (y[i] - sum) / U(i, i);
        // x[i] = (y[i] - sum) / U(i, i);
    }

    return x;
}

std::vector<std::complex<double>> LUSolveLinearSystem(
    const std::tuple<Matrix<std::complex<double>>,
                     Matrix<std::complex<double>>,
                     std::vector<int>>& l_u_p,
    const std::vector<std::complex<double>>& b) {
    // Matrix<std::complex<double>> L;
    // Matrix<std::complex<double>> U;
    // std::vector<int> perm;

    auto [L, U, perm] = l_u_p;

    auto bv = b;
    int n = b.size();

    // Perform LU decomposition (manual implementation)

    for (int i = 0; i < n; ++i) { bv[i] = b[perm[i]]; }

    // Forward substitution (solve Ly = b)
    std::vector<std::complex<double>> y(n);
    for (int i = 0; i < n; ++i) {
        std::complex<double> sum = 0.0;
        for (int j = 0; j < i; ++j) { sum += L(i, j) * y[j]; }
        //       for (int j = 0; j < i; ++j) { sum += L(i, j) * y[j]; }
        y[i] = bv[i] - sum;
    }

    // Backward substitution (solve Ux = y)
    std::vector<std::complex<double>> x(n);
    for (int i = n - 1; i >= 0; --i) {
        std::complex<double> sum = 0.0;
        for (int j = i + 1; j < n; ++j) { sum += U(i, j) * x[j]; }
        x[i] = (y[i] - sum) / U(i, i);
        // x[i] = (y[i] - sum) / U(i, i);
    }

    return x;
}

// // Newton trace iteration for NLEP
// std::pair<std::complex<double>, Matrix<std::complex<double>>>
// NewtonTraceIteration(std::complex<double> lambda, double tol) {
//     int max_iter = 100;
//     for (int i = 0; i < max_iter; ++i) {
//         Matrix<std::complex<double>> F_lambda = F(lambda);
//         Matrix<std::complex<double>>
//         linear_solution_matrix(F_lambda.getRows(),
//                                                             F_lambda.getCols());

//         Matrix<std::complex<double>> J_lambda = Jacobian(lambda);

//         auto l_u_p = F_lambda.luDecomposition();

//         for (int j = 0; j < F_lambda.getCols(); ++j) {
//             linear_solution_matrix.setCol(
//                 j, LUSolveLinearSystem(l_u_p, J_lambda.getCol(j)));
//         };  // getCol and setCol is slow, need to be optimized;

//         // Update eigenvalue and eigenvector
//         auto one_over_linear_solution_matrix_trace =
//             1.0 / linear_solution_matrix.trace();

//         if (std::abs(one_over_linear_solution_matrix_trace) >
//             std::abs(tol * lambda)) {
//             lambda -= one_over_linear_solution_matrix_trace;
//         } else {
//             return {lambda, F(lambda)};
//         }
//     }

//     // // Throw an exception or return a flag if convergence is not achieved
//     // throw std::runtime_error(
//     //     "Newton trace iteration did not converge within maximum
//     iterations");

//     return {lambda, F(lambda)};
// }

// Newton trace iteration for NLEP
std::pair<std::complex<double>, Matrix<std::complex<double>>>
NewtonTraceIterationSecantMethod(std::complex<double> lambda,
                                 const double& tol,
                                 Parameters& para,
                                 const Matrix<double>& coeff_matrix,
                                 const Grid<double>& grid_info,
                                 const int& max_iter) {
    Matrix<std::complex<double>> F_lambda = F(
        para.tau, lambda,
        [&para](double eta, double eta_p, std::complex<double> omega) {
            return para.kappa_f_tau(eta, eta_p, omega);
        },
        coeff_matrix, grid_info);

    Matrix<std::complex<double>> F_old_lambda = F_lambda;
    Matrix<std::complex<double>> linear_solution_matrix(F_lambda.getRows(),
                                                        F_lambda.getCols());

    Matrix<std::complex<double>> J_lambda =
        (F_lambda -
         F(
             para.tau, lambda * 0.99,
             [&para](double eta, double eta_p, std::complex<double> omega) {
                 return para.kappa_f_tau(eta, eta_p, omega);
             },
             coeff_matrix, grid_info)) /
        (0.01 * lambda);

    for (int i = 0; i < max_iter; ++i) {
        auto l_u_p = F_lambda.luDecomposition();

        for (int j = 0; j < F_lambda.getCols(); ++j) {
            linear_solution_matrix.setCol(
                j, LUSolveLinearSystem(l_u_p, J_lambda.getCol(j)));
        };  // getCol and setCol is slow, need to be optimized;

        // Update eigenvalue and eigenvector
        auto d_lambda = 1.0 / linear_solution_matrix.trace();

        if (std::abs(d_lambda) > std::abs(tol * lambda)) {
            lambda -= d_lambda;
            std::cout << lambda << std::endl;
        } else {
            return {lambda, F_lambda};
        }

        F_old_lambda = F_lambda;
        F_lambda = F(
            para.tau, lambda,
            [&para](double eta, double eta_p, std::complex<double> omega) {
                return para.kappa_f_tau(eta, eta_p, omega);
            },
            coeff_matrix, grid_info);
        J_lambda = (F_old_lambda - F_lambda) / (d_lambda);
    }

    // // Throw an exception or return a flag if convergence is not achieved
    // throw std::runtime_error(
    //     "Newton trace iteration did not converge within maximum iterations");

    return {lambda, F_lambda};
}

// int main() {
//     // Define the initial guess for the eigenvector (lambda0)
//     std::complex<double> lambda0 = 2.4;

//     // Set the tolerance for convergence
//     double tol = 1e-6;

//     // try {
//     std::pair<std::complex<double>, Matrix<std::complex<double>>> result =
//         NewtonTraceIterationSecantMethod(lambda0, tol);

//     // Print the resulting eigenvalue and (normalized) eigenvector
//     // std::complex<double> norm = 0.0;
//     // for (const auto& val : result.second) { norm += std::norm(val); }

//     std::cout << "Eigenvalue: " << result.first << std::endl;

//     // for (const auto& val : result.second) { std::cout << val << " "; }
//     // std::cout << std::endl;
//     // Matrix<std::complex<double>> test(2, 2);
//     // test.setRow(0, {{1.0, 0.5}, {2.0, 0.1}});
//     // test.setRow(1, {{1.3, 0.2}, {2.0, 0.3}});

//     // for (const auto& val : result.second - test) { std::cout << val << "
//     "; }
//     // std::cout << std::endl;

//     // };

//     return 0;
// }
