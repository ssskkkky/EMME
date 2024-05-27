#include <complex>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <utility>

#include "Grid.h"
#include "JsonParser.h"
#include "Matrix.h"
#include "Parameters.h"
// #include "fenv.h" //this is for check inf or nan
#include "functions.h"
#include "singularity_handler.h"
#include "solver.h"

int main() {
    std::string filename = "input.json";

    auto input = util::json::parse_file(filename);

    std::complex<double> omega_initial_guess(input["initial_guess"][0],
                                             input["initial_guess"][1]);

    // reserve stack space
    alignas(Stellarator) std::byte buffer[sizeof(Stellarator)];

    Parameters* para_ptr = nullptr;
    // Parameters are Stellarator are both trivially destructible, no need to
    // bother calling their destructors.
    if (!std::string{"tokamak"}.compare(input["conf"])) {
        para_ptr = new (buffer) Parameters(
            input["q"], input["shat"], input["tau"], input["epsilon_n"],
            input["eta_i"], input["eta_e"], input["b_theta"], input["beta_e"],
            input["R"], input["vt"], input["length"], input["theta"],
            input["npoints"], input["iteration_step_limit"]);
    } else if (!std::string{"stellarator"}.compare(input["conf"])) {
        para_ptr = new (buffer) Stellarator(
            input["q"], input["shat"], input["tau"], input["epsilon_n"],
            input["eta_i"], input["eta_e"], input["b_theta"], input["beta_e"],
            input["R"], input["vt"], input["length"], input["theta"],
            input["npoints"], input["iteration_step_limit"], input["eta_k"],
            input["lh"], input["mh"], input["epsilon_h_t"], input["alpha_0"],
            input["r_over_R"]);
    } else {
        throw std::runtime_error("Input configuration not supported yet.");
    }
    auto& para = *para_ptr;

    auto length = para.length;
    auto npoints = para.npoints;

    Grid<double> grid_info(length, npoints);

    Matrix<double> coeff_matrix = SingularityHandler(npoints);

    std::ofstream outfile("emme_eigen_vector.csv");
    std::ofstream eigenvalue("emme_eigen_value.csv");

    double tol = 1e-6;

    for (unsigned int i = 0; i <= 0; i++) {
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
