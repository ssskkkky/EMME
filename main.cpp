#include <complex>
#include <fstream>
#include <iostream>
#include <utility>

#include "Grid.h"
#include "JsonParser.h"
#include "Matrix.h"
#include "Parameters.h"
#include "functions.h"
#include "singularity_handler.h"
#include "solver.h"

int main() {
    std::string filename = "input.json";

    auto input = util::json::parse_file(filename);

    std::complex<double> omega_initial_guess(input["initial_guess"][0],
                                             input["initial_guess"][1]);

    Parameters para(input["q"], input["shat"], input["tau"], input["epsilon_n"],
                    input["eta_i"], input["eta_e"], input["b_theta"],
                    input["beta_e"], input["R"], input["vt"], input["length"],
                    input["theta"], input["npoints"],
                    input["iteration_step_limit"]);

    auto length = para.length;
    auto npoints = para.npoints;

    Grid<double> grid_info(length, npoints);

    // Matrix iter_matrix =
    //     matrix_assembler(tau, omega_iter, kappa_f, coeff_matrix,
    //     grid_info);

    // Matrix coeff_matrix =
    //     singularity_handler(grid_info, gauss_order, interpolation_order);

    Matrix<double> coeff_matrix = SingularityHandler(npoints);

    // IterateSolver iter_solver(tau, omega_initial, matrix_assembler,
    //                           coeff_matrix, grid_info, iteration_step_limit);

    // iter_solver.run();

    // Vector eigen_vector = eigen_vector_solver(iter_solver.matrix);

    // output(filename);

    // Set the tolerance for convergence
    // std::pair<std::complex<double>, Matrix<std::complex<double>>> result;
    std::ofstream outfile("emme_eigen_vector.csv");
    std::ofstream eigenvalue("emme_eigen_value.csv");

    double tol = 1e-6;
    for (unsigned int i = 0; i <= 40; i++) {
        std::cout << para.q << std::endl;
        auto result = NewtonTraceIterationSecantMethod(
            omega_initial_guess, tol, para, coeff_matrix, grid_info,
            para.iteration_step_limit);

        std::cout << "Eigenvalue: " << result.first.real() << " "
                  << result.first.imag() << std::endl;
        eigenvalue << result.first.real() << " " << result.first.imag()
                   << std::endl;
        auto null_space = NullSpace(result.second, 1e-1);
        if (!outfile.is_open()) {
            // Handle error
            return 1;
        }

        outfile << null_space;
        para.q += 0.05;
        omega_initial_guess = result.first;
    }

    outfile.close();

    // std::cout << para.kappa_f_tau(3.0, 2.0, 1.0);

    // std::cout << 1.0 / para.lambda_f_tau(3.0, 2.0, 1.0) *
    //                  util::bessel_ic(
    //                      0, std::sqrt((para.bi(3.0) * para.bi(2.0)) /
    //                                   para.lambda_f_tau(3.0, 2.0, 1.0)))
    //                                   *
    //                  std::exp(-(para.bi(3.0) + para.bi(2.0)) /
    //                           (2.0 * para.lambda_f_tau(3.0, 2.0, 1.0)));

    // std::cout << 1.0 / para.lambda_f_tau(3.0, 2.0, 1.0);

    // std::cout << std::sqrt((para.bi(3.0) * para.bi(2.0))) /
    //                  para.lambda_f_tau(3.0, 2.0, 1.0);

    // std::cout << util::bessel_ic(0,
    //                              std::sqrt((para.bi(3.0) * para.bi(2.0))
    //                              /
    //                                        para.lambda_f_tau(3.0, 2.0, 1.0)));

    // std::cout << std::exp(-(para.bi(3.0) + para.bi(2.0)) /
    //                       (2.0 * para.lambda_f_tau(3.0, 2.0, 1.0)));

    // std::cout << para.bi(3.0) * para.bi(2.0);

    // std::cout << para.bi(3.0);
    // std::cout << para.bi(2.0);

    return 0;
}
