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

using namespace util::json;

int main() {
    Timer::get_timer().start_timing("All");
    using namespace std::string_literals;
    std::string filename = "input.json";
    auto input_all = util::json::parse_file(filename);
    auto input = util::json::parse_file(filename);
    std::ofstream outfile("emme_eigen_vector.csv",
                          input["output"].as_string() == std::string{"append"}
                              ? std::ios::app
                              : std::ios::trunc);
    std::ofstream eigenvalue("emme_eigen_value.csv",
                             input["output"].as_string() == "append"
                                 ? std::ios::app
                                 : std::ios::trunc);
    double tol = input["iteration_precision"];
    std::complex<double> omega_initial_guess(input["initial_guess"][0],
                                             input["initial_guess"][1]);

    // record omega_head as the initial guess for other side scanning.
    auto omega_head = omega_initial_guess;

    alignas(Stellarator) std::byte buffer[sizeof(Stellarator)];

    // create output object
    auto result = Value::create_object();
    result["input"] = input_all.clone();
    // find out which parameter is setup for scan
    std::string scan_key{};
    Value scan_opt;
    for (auto& [key, val] : input_all.as_object()) {
        if (val.is_object()) {
            scan_key = key;
            scan_opt = val.clone();
        }
    }

    Timer::get_timer().start_timing("initial");
    auto head = scan_opt["head"];
    auto tail_array = scan_opt["tail"].as_array();
    for (unsigned int ii = 0; ii < 2; ii++) {
        auto tail = tail_array[ii];
        if (ii > 0) {
            head = scan_opt["head"] +
                   std::copysign(scan_opt["step"].as_number<double>(),
                                 (tail - scan_opt["head"]));
            omega_initial_guess = omega_head;
        }
        for (input[scan_key] = head.as_number<double>();
             std::abs(input[scan_key].as_number<double>() - scan_opt["head"]) <=
             (std::abs(tail - scan_opt["head"]) +
              // This stupid but worked 0.01 term is for dealing with
              // the float error
              0.01 * std::abs(scan_opt["step"].as_number<double>()));
             input[scan_key] += std::copysign(
                 scan_opt["step"].as_number<double>(), (tail - head))) {
            Parameters* para_ptr = nullptr;
            // Parameters and Stellarator are both trivially
            // destructible, no need to bother calling their
            // destructors.
            if (!std::string{"tokamak"}.compare(input["conf"])) {
                para_ptr = new (buffer) Parameters(
                    input["q"], input["shat"], input["tau"], input["epsilon_n"],
                    input["eta_i"], input["eta_e"], input["b_theta"],
                    input["beta_e"], input["R"], input["vt"], input["length"],
                    input["theta"], input["npoints"],
                    input["iteration_step_limit"],
                    input["integration_precision"],
                    input["integration_accuracy"],
                    input["integration_iteration_limit"],
                    input["integration_start_points"], input["arc_coeff"]);
            } else if (!std::string{"stellarator"}.compare(input["conf"])) {
                para_ptr = new (buffer) Stellarator(
                    input["q"], input["shat"], input["tau"], input["epsilon_n"],
                    input["eta_i"], input["eta_e"], input["b_theta"],
                    input["beta_e"], input["R"], input["vt"], input["length"],
                    input["theta"], input["npoints"],
                    input["iteration_step_limit"],
                    input["integration_precision"],
                    input["integration_accuracy"],
                    input["integration_iteration_limit"],
                    input["integration_start_points"], input["arc_coeff"],
                    input["eta_k"], input["lh"], input["mh"],
                    input["epsilon_h_t"], input["alpha_0"], input["r_over_R"]);
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
            std::cout << scan_key << ":" << input[scan_key] << std::endl;
            Timer::get_timer().pause_timing("initial");

            for (int j = 0; j <= para.iteration_step_limit; j++) {
                Timer::get_timer().start_timing("newtonTracSecantIteration");
                eigen_solver.newtonTraceSecantIteration();
                Timer::get_timer().pause_timing("newtonTracSecantIteration");

                std::cout << eigen_solver.eigen_value << std::endl;
                if (std::abs(eigen_solver.d_eigen_value) <
                    std::abs(tol * eigen_solver.eigen_value)) {
                    break;
                }
            }

            std::cout << "Eigenvalue: " << eigen_solver.eigen_value.real()
                      << " " << eigen_solver.eigen_value.imag() << std::endl;
            eigenvalue << eigen_solver.eigen_value.real() << " "
                       << eigen_solver.eigen_value.imag() << std::endl;
            Timer::get_timer().start_timing("SVD");
            auto null_space = eigen_solver.nullSpace();
            Timer::get_timer().pause_timing("SVD");
            if (!outfile.is_open()) {
                // Handle error
                return 1;
            }
            outfile << null_space << std::endl;
            flush(eigenvalue);
            flush(outfile);
            omega_initial_guess = eigen_solver.eigen_value;
            if (ii == 0 && input[scan_key] == head) {
                omega_head = omega_initial_guess;
            }
        }
    }

    outfile.close();
    Timer::get_timer().pause_timing("All");
    Timer::get_timer().print();
    std::cout << std::endl;

    return 0;
}
