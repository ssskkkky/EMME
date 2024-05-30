#include <complex>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <utility>

#include "Grid.h"
#include "JsonParser.h"
#include "Matrix.h"
#include "Parameters.h"
#include "Timer.h"
#include "functions.h"
#include "singularity_handler.h"
#include "solver.h"

int main() {
    Timer::get_timer().start_timing("All");

    using namespace std::string_literals;
    std::string filename = "input.json";
    auto input_all = util::json::parse_file(filename);
    auto output_mode = input_all["output"].as_string() == std::string{"append"}
                           ? std::ios::app
                           : std::ios::trunc;
    std::ofstream outfile("emme_eigen_vector.csv", output_mode);
    std::ofstream eigenvalue("emme_eigen_value.csv", output_mode);
    double tol = input_all["iteration_precision"];

    auto input = input_all.clone();
    std::complex<double> omega_initial_guess(input["initial_guess"][0],
                                             input["initial_guess"][1]);
    // initial_guess is not designed for scanning.

    // reserve stack space for parameter
    alignas(Stellarator) std::byte buffer[sizeof(Stellarator)];

    std::string scan_key;
    util::json::Value scan_opt;

    // check which parameter is set for scan
    for (auto& [key, val] : input_all.as_object()) {
        if (val.is_object()) {
            scan_key = key;
            scan_opt = val.clone();
            break;
        }
    }

    for (input[scan_key] = scan_opt["head"];
         std::abs(input[scan_key] - scan_opt["head"]) <=
         (std::abs(scan_opt["tail"] - scan_opt["head"]) +
          // This stupid but worked 0.01 term is for dealing with the
          // float error
          0.01 * std::abs(scan_opt["step"]));
         input[scan_key] += std::copysign(
             scan_opt["step"], (scan_opt["tail"] - scan_opt["head"]))) {
        Timer::get_timer().start_timing("initial");
        Parameters* para_ptr = nullptr;
        // Parameters and Stellarator are both trivially destructible,
        // no need to bother calling their destructors.
        if (!std::string{"tokamak"}.compare(input["conf"])) {
            para_ptr = new (buffer) Parameters(
                input["q"], input["shat"], input["tau"], input["epsilon_n"],
                input["eta_i"], input["eta_e"], input["b_theta"],
                input["beta_e"], input["R"], input["vt"], input["length"],
                input["theta"], input["npoints"], input["iteration_step_limit"],
                input["integration_precision"], input["integration_accuracy"],
                input["integration_iteration_limit"],
                input["integration_start_points"]);
        } else if (!std::string{"stellarator"}.compare(input["conf"])) {
            para_ptr = new (buffer) Stellarator(
                input["q"], input["shat"], input["tau"], input["epsilon_n"],
                input["eta_i"], input["eta_e"], input["b_theta"],
                input["beta_e"], input["R"], input["vt"], input["length"],
                input["theta"], input["npoints"], input["iteration_step_limit"],
                input["integration_precision"], input["integration_accuracy"],
                input["integration_iteration_limit"],
                input["integration_start_points"], input["eta_k"], input["lh"],
                input["mh"], input["epsilon_h_t"], input["alpha_0"],
                input["r_over_R"]);
        } else {
            throw std::runtime_error("Input configuration not supported yet.");
        }
        auto& para = *para_ptr;

        auto length = para.length;
        auto npoints = para.npoints;

        Grid<double> grid_info(length, npoints);

        Matrix<double> coeff_matrix = SingularityHandler(npoints);

        auto eigen_solver = EigenSolver<Matrix<std::complex<double>>>(
            para, omega_initial_guess, coeff_matrix, grid_info);
        std::cout << scan_key << ":" << input[scan_key] << std::endl;
        Timer::get_timer().pause_timing("initial");

        for (int j = 0; j <= para.iteration_step_limit; j++) {
            Timer::get_timer().start_timing("newtonTracSecantIteration");
            eigen_solver.newtonTraceSecantIteration();
            Timer::get_timer().pause_timing("newtonTracSecantIteration");

            std::cout << eigen_solver.eigen_value << '\n';
            if (std::abs(eigen_solver.d_eigen_value) <
                std::abs(tol * eigen_solver.eigen_value)) {
                break;
            }
        }

        std::cout << "Eigenvalue: " << eigen_solver.eigen_value.real() << " "
                  << eigen_solver.eigen_value.imag() << '\n';
        eigenvalue << eigen_solver.eigen_value << '\n';

        Timer::get_timer().start_timing("SVD");
        auto null_space = eigen_solver.nullSpace();
        Timer::get_timer().pause_timing("SVD");
        if (!outfile.is_open()) {
            std::cout << "Output file can not be opened.\n";
            return 1;
        }
        outfile << null_space;
        omega_initial_guess = eigen_solver.eigen_value;
    }

    outfile.close();
    Timer::get_timer().pause_timing("All");
    Timer::get_timer().print();

    std::cout << '\n';

    return 0;
}
