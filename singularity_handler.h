#ifndef SINGULARITYHANDLER_H
#define SINGULARITYHANDLER_H

#include <vector>
#include "Matrix.h"

Matrix<double> SingularityHandler(int n) {
    std::vector<double> coeff = {0.0,
                                 2.951388888888883,
                                 -2.4305555555555305,
                                 4.166666666667441,
                                 -0.3472222222224549,
                                 1.159722222222284};
    Matrix<double> coeff_matrix(n, n);
    for (int i; i < n; i++) {
        for (int j; j < n; j++) {
            int diff = std::abs(i - j);
            if (diff <= 5) {
                coeff_matrix(i, j) = coeff[diff];
            } else {
                coeff_matrix(i, j) = 1.0;
            }
            if (j == 0 || j == n - 1) { coeff_matrix(i, j) -= 0.5; }
        }
    }

    return coeff_matrix;
}

#endif
