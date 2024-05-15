#ifndef SOLVER_H  // Replace MATRIX_H with your unique guard macro name
#define SOLVER_H

#include <iostream>

#include <complex>
#include <vector>

#include "Matrix.h"

std::pair<std::complex<double>, Matrix<std::complex<double>>>
NewtonTraceIteration(std::complex<double> lambda, double tol);

std::pair<std::complex<double>, Matrix<std::complex<double>>>
NewtonTraceIterationSecantMethod(std::complex<double> lambda, double tol);

template <typename T>
Matrix<T> nullSpace(const Matrix<T>& A) {
    // Check if the matrix is square
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Input matrix must be square.");
    }

    // Perform Singular Value Decomposition (SVD)
    // You'll need an external library like Eigen or LAPACK for SVD
    // Replace this placeholder with your preferred SVD implementation
    // which returns the singular values (S), left singular vectors (U),
    // and right singular vectors (V)
    std::vector<T> S;
    Matrix<T> U;
    Matrix<T> V;
    // ... perform SVD on A and store results in S, U, V

    // Identify null space basis vectors and construct the null space matrix
    int null_space_dim = 0;  // Count the number of null space basis vectors
    for (int i = 0; i < A.cols(); ++i) {
        // Check for singular values close to zero (tolerance approach)
        if (std::abs(S[i]) < std::numeric_limits<T>::epsilon()) {
            null_space_dim++;
        }
    }

    // Create a result matrix to store the null space basis vectors as columns
    Matrix<T> null_space(A.rows(), null_space_dim);
    int col_index = 0;
    for (int i = 0; i < A.cols(); ++i) {
        // Check for singular values close to zero (tolerance approach)
        if (std::abs(S[i]) < std::numeric_limits<T>::epsilon()) {
            // Extract the corresponding right singular vector and store in null
            // space matrix
            null_space.setCol(col_index, U.getCol(i));
            col_index++;
        }
    }

    return null_space;
};

#endif
