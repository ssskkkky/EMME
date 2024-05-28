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
    std::ofstream outfile("emme_eigen_vector.csv");
    std::ofstream eigenvalue("emme_eigen_value.csv");
    double tol = 1e-6;

    using namespace std::string_literals;
    auto input_all = util::json::parse_file(filename);
    auto input = util::json::parse_file(filename);
    std::complex<double> omega_initial_guess(
        input["initial_guess"s][0],
        input["initial_guess"s]
             [1]);  // initial_guess is not designed for scanning.
    // reserve stack space
    alignas(Stellarator) std::byte buffer[sizeof(Stellarator)];

    for (auto& [key, val] : input_all.as_object()) {
        if (val.is_object()) {
            for (input[key] = val["head"s].as_number<double>();
                 input[key] <= val["tail"s];
                 input[key] += val["step"s].as_number<double>()) {
                Parameters* para_ptr = nullptr;
                // Parameters and Stellarator are both trivially destructible,
                // no need to bother calling their destructors.
                if (!std::string{"tokamak"}.compare(input["conf"s])) {
                    para_ptr = new (buffer) Parameters(
                        input["q"s], input["shat"s], input["tau"s],
                        input["epsilon_n"s], input["eta_i"s], input["eta_e"s],
                        input["b_theta"s], input["beta_e"s], input["R"s],
                        input["vt"s], input["length"s], input["theta"s],
                        input["npoints"s], input["iteration_step_limit"s]);
                } else if (!std::string{"stellarator"}.compare(
                               input["conf"s])) {
                    para_ptr = new (buffer) Stellarator(
                        input["q"s], input["shat"s], input["tau"s],
                        input["epsilon_n"s], input["eta_i"s], input["eta_e"s],
                        input["b_theta"s], input["beta_e"s], input["R"s],
                        input["vt"s], input["length"s], input["theta"s],
                        input["npoints"s], input["iteration_step_limit"s],
                        input["eta_k"s], input["lh"s], input["mh"s],
                        input["epsilon_h_t"s], input["alpha_0"s],
                        input["r_over_R"s]);
                } else {
                    throw std::runtime_error(
                        "Input configuration not supported yet.");
                }
                auto& para = *para_ptr;

                auto length = para.length;
                auto npoints = para.npoints;

                Grid<double> grid_info(length, npoints);

                Matrix<double> coeff_matrix = SingularityHandler(npoints);

                auto eigen_solver = EigenSolver<Matrix<std::complex<double>>>(
                    para, omega_initial_guess, coeff_matrix, grid_info);
                std::cout << key << ":" << input[key] << std::endl;

                for (int j = 0; j <= para.iteration_step_limit; j++) {
                    eigen_solver.newtonTraceSecantIteration();
                    std::cout << eigen_solver.eigen_value << std::endl;
                    if (std::abs(eigen_solver.d_eigen_value) <
                        std::abs(tol * eigen_solver.eigen_value)) {
                        break;
                    }
                }

                std::cout << "Eigenvalue: " << eigen_solver.eigen_value.real()
                          << " " << eigen_solver.eigen_value.imag()
                          << std::endl;
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
                // para.beta_e += 0.0005;
                // para.parameterInit();
                omega_initial_guess = eigen_solver.eigen_value;
            }
        }
    }
    outfile.close();

    return 0;
}
