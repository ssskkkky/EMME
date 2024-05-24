#include <complex>
#include <fstream>
#include <iostream>
#include <utility>

#include "Grid.h"
#include "Matrix.h"
#include "Parameters.h"
// #include "fenv.h" //this is for check inf or nan
#include "functions.h"
#include "singularity_handler.h"
#include "solver.h"

int main() {
    // feenableexcept(FE_DIVBYZERO);
    std::string filename = "emme.in";
    std::ifstream input_file(filename);

    double q_input;
    double shat_input;
    double tau_input;
    double epsilon_n_input;
    double eta_i_input;
    double eta_e_input;
    double b_theta_input;
    double beta_e_input;
    double R_input;
    double vt_input;
    double length_input;
    double theta_input;
    int npoints_input;
    int iteration_step_limit_input;
    double initial_guess_real;
    double initial_guess_imag;

    input_file >> q_input >> shat_input >> tau_input >> epsilon_n_input >>
        eta_i_input >> eta_e_input >> b_theta_input >> beta_e_input >>
        R_input >> vt_input >> length_input >> theta_input >> npoints_input >>
        iteration_step_limit_input >> initial_guess_real >> initial_guess_imag;

    input_file.close();

    // const auto [q; shat, tau, epsilon_n, eta_i, b_theta, R, vt,
    // length, theta,
    //             npoints, iteration_step_limit] =
    //             readInputFromFile(filename);

    // std::tie(omega_s_i, omega_d_bar) = initiallize(vt, tau);

    // Parameters para(q_input, shat_input, tau_input, epsilon_n_input,
    //                 eta_i_input, b_theta_input, R_input, vt_input,
    //                 length_input, theta_input, npoints_input,
    //                 iteration_step_limit_input);

    std::complex<double> omega_initial_guess(initial_guess_real,
                                             initial_guess_imag);

    Parameters para(q_input, shat_input, tau_input, epsilon_n_input,
                    eta_i_input, eta_e_input, b_theta_input, beta_e_input,
                    R_input, vt_input, length_input, theta_input, npoints_input,
                    iteration_step_limit_input);

    auto length = para.length;
    auto npoints = para.npoints;

    Grid<double> grid_info(length, npoints);

    Matrix<double> coeff_matrix = SingularityHandler(npoints);

    std::ofstream outfile("emme_eigen_vector.csv");
    std::ofstream eigenvalue("emme_eigen_value.csv");

    double tol = 1e-6;

    for (unsigned int i = 0; i <= 40; i++) {
        auto eigen_solver = EigenSolver<Matrix<std::complex<double>>>(
            para, omega_initial_guess, coeff_matrix, grid_info);
        std::cout << eigen_solver.para.q << std::endl;

        for (int j = 0; j <= para.iteration_step_limit; j++) {
            eigen_solver.newtonTraceSecantIteration();
            std::cout << eigen_solver.eigen_value << std::endl;
            if (std::abs(eigen_solver.d_eigen_value) <
                std::abs(tol * eigen_solver.eigen_value)) {
                break;
            }
        }

        std::cout << "Eigenvalue: " << eigen_solver.eigen_value.real() << " "
                  << eigen_solver.eigen_value.imag() << std::endl;
        eigenvalue << eigen_solver.eigen_value.real() << " "
                   << eigen_solver.eigen_value.imag() << std::endl;
        auto null_space = eigen_solver.nullSpace();
        if (!outfile.is_open()) {
            // Handle error
            return 1;
        }
        outfile << null_space;
        flush(eigenvalue);
        flush(outfile);
        para.q += 0.05;
        para.parameterInit();
        omega_initial_guess = eigen_solver.eigen_value;
    }

    outfile.close();

    return 0;
}
