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

using namespace util::json;

auto solve_once(auto& input,
                auto& result_array,
                const auto& scan_key,
                auto omega_initial_guess) {
    auto& timer = Timer::get_timer();
    double tol = input["iteration_precision"];

    timer.start_timing("initial");
    alignas(Stellarator) std::byte buffer[sizeof(Stellarator)];
    Parameters* para_ptr = nullptr;
    // Parameters and Stellarator are both trivially
    // destructible, no need to bother calling their
    // destructors.
    // TODO: use factory method pattern
    if (!std::string{"tokamak"}.compare(input["conf"])) {
        para_ptr = new (buffer) Parameters(input);
    } else if (!std::string{"stellarator"}.compare(input["conf"])) {
        para_ptr = new (buffer) Stellarator(input);
    } else if (!std::string{"cylinder"}.compare(input["conf"])) {
        para_ptr = new (buffer) Cylinder(input);
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
    timer.pause_timing("initial");

    for (int j = 0; j <= para.iteration_step_limit; j++) {
        timer.start_timing("newtonTracSecantIteration");
        eigen_solver.newtonTraceSecantIteration();
        timer.pause_timing("newtonTracSecantIteration");

        std::cout << eigen_solver.eigen_value << std::endl;
        if (std::abs(eigen_solver.d_eigen_value) <
            std::abs(tol * eigen_solver.eigen_value)) {
            break;
        }
    }

    std::cout << "Eigenvalue: " << eigen_solver.eigen_value.real() << " "
              << eigen_solver.eigen_value.imag() << '\n';

    // store eigenvalue and eigenvector to result
    auto result_unit = Value::create_object();
    result_unit["scan_value"] = input[scan_key];
    auto& eva = result_unit["eigenvalue"] = Value::create_array(2);
    eva[0] = eigen_solver.eigen_value.real();
    eva[1] = eigen_solver.eigen_value.imag();

    timer.start_timing("SVD");
    result_unit["eigenvector"] =
        Value::create_typed_array(eigen_solver.nullSpace());
    timer.pause_timing("SVD");

    result_array.push_back(std::move(result_unit));

    return eigen_solver.eigen_value;
}

int main() {
    auto& timer = Timer::get_timer();
    timer.start_timing("All");
    using namespace std::string_literals;

    std::string filename = "input.json";
    std::string output_filename = "output.json";

    auto input_all = util::json::parse_file(filename);
    auto input = input_all.clone();
    std::complex<double> omega_initial_guess(input["initial_guess"][0],
                                             input["initial_guess"][1]);

    // record omega_head as the initial guess for other side scanning.
    auto omega_head = omega_initial_guess;

    // create output object
    auto result = Value::create_object();
    result["input"] = input_all.clone();

    // build info from preprocessor invocation command arguments, works when
    // using GNU Makefiles only

#ifdef EMME_COMMIT_HASH
    result["git_commit_hash"] = EMME_COMMIT_HASH;
#endif
#ifdef EMME_BUILD_DATE
    result["build_time"] = EMME_BUILD_DATE;
#endif
    result["run_time"] = util::get_date_string();

    // find out which parameter is setup for scan
    std::string scan_key{};
    Value scan_opt;
    for (auto& [key, val] : input_all.as_object()) {
        if (val.is_object()) {
            scan_key = key;
            scan_opt = val.clone();
        }
    }

    result["scan_parameter"] = scan_key;
    result["result"] = Value::create_array();
    auto& result_array = result["result"].as_array();

    auto head = scan_opt["head"];
    const auto step = scan_opt["step"];
    auto tail_array = scan_opt["tail"].as_array();
    for (unsigned int ii = 0; ii < 2; ii++) {
        const auto& tail = tail_array[ii];
        for (input[scan_key] = head;
             std::abs(input[scan_key] - scan_opt["head"]) <=
             (std::abs(tail - scan_opt["head"]) +
              // This stupid but worked 0.01 term is for dealing with
              // the float error
              0.01 * std::abs(step));
             input[scan_key] += std::copysign(step, (tail - head))) {
            std::cout << scan_key << ":" << input[scan_key] << '\n';
            omega_initial_guess =
                solve_once(input, result_array, scan_key, omega_initial_guess);

            // store the result of head
            if (ii == 0 && input[scan_key] == head) {
                omega_head = omega_initial_guess;
            }
        }

        // After scanning one direction, reset the head to middle, but one step
        // towards another direction, preventing recalculation.
        head = scan_opt["head"] +
               std::copysign(step, tail_array[1] - scan_opt["head"]);
        omega_initial_guess = omega_head;
    }
    timer.start_timing("Output");
    std::ofstream output(output_filename);
    output << result.dump();
    timer.pause_timing("Output");

    timer.pause_timing("All");
    timer.print();
    std::cout << '\n';

    return 0;
}
