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
#include "Timer.h"
#include "functions.h"
#include "singularity_handler.h"
#include "solver.h"

int main() {
    Timer::get_timer().start_timing("All");
    std::string filename = "input.json";
    std::ofstream outfile("emme_eigen_vector.csv");
    std::ofstream eigenvalue("emme_eigen_value.csv");
    double tol = 1e-6;

    using namespace std::string_literals;
    auto input_all = util::json::parse_file(filename);
    auto input = util::json::parse_file(filename);
    std::complex<double> omega_initial_guess(
        input["initial_guess"][0],
        input["initial_guess"]
             [1]);  // initial_guess is not designed for scanning.
    // reserve stack space
    alignas(Stellarator) std::byte buffer[sizeof(Stellarator)];

    for (auto& [key, val] : input_all.as_object()) {
        Timer::get_timer().start_timing("initial");
        if (val.is_object()) {
            for (input[key] = val["head"].as_number<double>();
                 std::abs(input[key].as_number<double>() - val["head"]) <=
                 std::abs(val["tail"].as_number<double>() - val["head"]);
                 input[key] += val["step"].as_number<double>()) {
                Parameters* para_ptr = nullptr;
                // Parameters and Stellarator are both trivially destructible,
                // no need to bother calling their destructors.
                if (!std::string{"tokamak"}.compare(input["conf"])) {
                    para_ptr = new (buffer) Parameters(
                        input["q"], input["shat"], input["tau"],
                        input["epsilon_n"], input["eta_i"], input["eta_e"],
                        input["b_theta"], input["beta_e"], input["R"],
                        input["vt"], input["length"], input["theta"],
                        input["npoints"], input["iteration_step_limit"],
                        input["integration_precision"],
                        input["integration_accuracy"],
                        input["integration_iteration_limit"],
                        input["integration_start_points"]);
                } else if (!std::string{"stellarator"}.compare(input["conf"])) {
                    para_ptr = new (buffer) Stellarator(
                        input["q"], input["shat"], input["tau"],
                        input["epsilon_n"], input["eta_i"], input["eta_e"],
                        input["b_theta"], input["beta_e"], input["R"],
                        input["vt"], input["length"], input["theta"],
                        input["npoints"], input["iteration_step_limit"],
                        input["integration_precision"],
                        input["integration_accuracy"],
                        input["integration_iteration_limit"],
                        input["integration_start_points"], input["eta_k"],
                        input["lh"], input["mh"], input["epsilon_h_t"],
                        input["alpha_0"], input["r_over_R"]);
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
                Timer::get_timer().pause_timing("initial");

                for (int j = 0; j <= para.iteration_step_limit; j++) {
                    Timer::get_timer().start_timing(
                        "newtonTracSecantIteration");
                    eigen_solver.newtonTraceSecantIteration();
                    Timer::get_timer().pause_timing(
                        "newtonTracSecantIteration");

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
                Timer::get_timer().start_timing("SVD");
                auto null_space = eigen_solver.nullSpace();
                Timer::get_timer().pause_timing("SVD");
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
    Timer::get_timer().pause_timing("All");
    Timer::get_timer().print();
    std::cout << std::endl;

    return 0;
}
