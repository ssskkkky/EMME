#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include "solver.h"

#include <complex>
#include <iostream>
#include <vector>

#include "Grid.h"
#include "Matrix.h"
#include "Parameters.h"
#include "functions.h"
#include "lapack.h"

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
