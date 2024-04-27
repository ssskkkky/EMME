#include <complex>
#include <iostream>

int main() {
    double omega_s_i, omega_d_bar;
    const auto [q, shat, tau, epsilon_n, eta_i, b_theta, R, vt, length, theta,
                npoints, iteration_step_limit] = read_input(filename);

    std::tie(omega_s_i, omega_d_bar) = initiallize(vt, tau);

    std::complex<double> omega_initial_guess(1.0, 1.0);

    Grid grid_info(length, npoints, theta);

    // Matrix iter_matrix =
    //     matrix_assembler(tau, omega_iter, kappa_f, coeff_matrix, grid_info);

    Matrix coeff_matrix =
        singularity_handler(grid_info, gauss_order, interpolation_order);

    IterateSolver iter_solver(tau, omega_initial, matrix_assembler,
                              coeff_matrix, grid_info, iteration_step_limit);

    iter_solver.run();

    Vector eigen_vector = eigen_vector_solver(iter_solver.matrix);

    output(filename);

    return 0;
}
