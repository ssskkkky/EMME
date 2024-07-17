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
#include "solver_pic.h"

using namespace util::json;

auto solve_once(auto& input,
                auto& omega_initial_guess,
                std::ofstream& eigen_matrix_file) {
    auto& timer = Timer::get_timer();
    double tol = input["iteration_precision"];

    timer.start_timing("initial");

    auto& para = Parameters::generate(input);

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

        std::cout << "        " << eigen_solver.eigen_value << '\n';
        if (std::abs(eigen_solver.d_eigen_value) <
            std::abs(tol * eigen_solver.eigen_value)) {
            break;
        }
    }

    std::cout << "        Eigenvalue: " << eigen_solver.eigen_value << '\n';
    timer.start_timing("Output");
    auto& v_output = eigen_solver.eigen_matrix;
    eigen_matrix_file.write(reinterpret_cast<char*>(v_output.data()),
                            sizeof(v_output(0, 0)) * v_output.size());

    // store eigenvalue and eigenvector to result

    auto single_result = Value::create_object();
    auto& eva = single_result["eigenvalue"] = Value::create_array(2);
    eva[0] = eigen_solver.eigen_value.real();
    eva[1] = eigen_solver.eigen_value.imag();
    timer.pause_timing("Output");

    timer.start_timing("SVD");
    single_result["eigenvector"] =
        Value::create_typed_array(eigen_solver.nullSpace());
    timer.pause_timing("SVD");

    omega_initial_guess = eigen_solver.eigen_value;
    return single_result;
}

auto solve_once_pic(const auto& input) {
    auto& para = Parameters::generate(input);
    const std::size_t marker_per_cell = 1 << 12;
    PIC_State<double> state(para, marker_per_cell);
    Integrator integrator(state);

    const std::size_t time_step = 160;
    const double dt = 0.25;

    Matrix<std::complex<double>> phi_history(time_step, para.npoints);
    for (std::size_t idx = 0; idx < time_step; ++idx) {
        integrator.step(dt);
        phi_history.setRow(idx, state.current_field());
        std::cout << idx + 1 << '/' << time_step << '\n';
    }

    std::ofstream ofs("PIC_phi.bin");
    ofs.write(reinterpret_cast<char*>(phi_history.data()),
              sizeof(phi_history(0, 0)) * phi_history.size());
}

auto get_scan_generator(const auto& para) {
    static double head, step, right_tail, left_tail, current, current_tail;
    static bool to_left, is_first;

    to_left = is_first = true;
    std::tie(head, step, left_tail, right_tail) = util::unpack(para);
    current = head;
    current_tail = left_tail;

    return [&]() mutable {
        // This stupid but worked 0.01 term is for dealing with the float
        // error
        auto within_range = [&]() {
            return std::abs(current - head) <=
                   (std::abs(current_tail - head) + 0.01 * std::abs(step));
        };

        if (!is_first) {
            current += std::copysign(step, (current_tail - head));
        }
        is_first = false;

        if (within_range()) {
            return std::make_tuple(true, false, current);
        } else {
            // go to another direction
            to_left = !to_left;
            current_tail = right_tail;
            current = head + std::copysign(step, (current_tail - head));

            return std::make_tuple(!to_left && within_range(), true, current);
        }
    };
}

auto filter_input(const auto& input_all) {
    auto input = input_all.clone();
    for (auto& [key, val] : input.as_object()) {
        if (val.is_object()) { val = val["head"]; }
    }
    return input;
}

int main() {
    std::string filename = "input.json";
    auto input = util::json::parse_file(filename);
    solve_once_pic(filter_input(input));
#if 0
    auto& timer = Timer::get_timer();
    timer.start_timing("All");
    using namespace std::string_literals;

    std::string filename = "input.json";
    std::string output_filename = "output.json";

    auto input_all = util::json::parse_file(filename);
    auto input = input_all.clone();
    std::complex<double> omega_initial_guess(input["initial_guess"][0],
                                             input["initial_guess"][1]);

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

    // scan_config["key"] = {heaad, step, tail, another_tail};
    // another_tail is optional, depending on the input
    std::unordered_map<std::string, std::array<double, 4>> scan_config;
    for (auto& [key, val] : input_all.as_object()) {
        if (val.is_object()) {
            auto iter = scan_config
                            .emplace(key, std::array<double, 4>{val["head"],
                                                                val["step"]})
                            .first;
            if (val["tail"].is_array()) {
                iter->second[2] = val["tail"][0];
                iter->second[3] = val["tail"][1];
            } else {
                iter->second[2] = val["tail"];
                iter->second[3] =
                    val["head"] +
                    .5 * std::copysign(val["step"], val["head"] - val["tail"]);
            }
        }
    }

    result["result"] = Value::create_object();
    auto& result_object = result["result"].as_object();

    if (scan_config.empty()) {
        // Do not need to scan any parameter
        auto result_unit = Value::create_object();
        result_unit["scan_key"] = "(None)";
        auto& scan_result_array = result_unit["scan_result"] =
            Value::create_array();
        std::cout << '\n';

        auto eigen_matrix_file_name = "eigenMatrics/eigenMatrix.bin";
        std::ofstream eigen_matrix_file(eigen_matrix_file_name,
                                        std::ios::binary);

        scan_result_array.as_array().push_back(
            solve_once(input_all, omega_initial_guess, eigen_matrix_file));
        result_object["(None)"] = std::move(result_unit);
    } else {
        auto omega = omega_initial_guess;
        for (const auto& [key, scan_para] : scan_config) {
            // for each parameter need to be scanned
            auto input = filter_input(input_all);
            auto get_scan_val = get_scan_generator(scan_para);
            auto [cont, turning, scan_value] = get_scan_val();
            auto result_unit = Value::create_object();
            result_unit["scan_key"] = key;
            auto& scan_value_array = result_unit["scan_values"] =
                Value::create_array();

            auto& scan_result_array = result_unit["scan_result"] =
                Value::create_array();
            std::cout << "\nScanning " << key << '\n';
            while (cont) {
                input[key] = scan_value;
                scan_value_array.as_array().push_back(scan_value);

                // begin to scan another direction
                if (turning) {
                    auto& new_omega = scan_result_array[0]["eigenvalue"];
                    omega.real(new_omega[0]);
                    omega.imag(new_omega[1]);
                }

                std::cout << "    " << key << ":" << scan_value << '\n';

                auto eigen_matrix_file_name = "eigenMatrics/" + key + "Eq" +
                                              std::to_string(scan_value) +
                                              ".bin";
                std::ofstream eigen_matrix_file(eigen_matrix_file_name,
                                                std::ios::binary);
                auto single_result =
                    solve_once(input, omega, eigen_matrix_file);
                single_result["eigenMatrix"] = eigen_matrix_file_name;
                single_result["scan_value"] = scan_value;
                scan_result_array.as_array().push_back(
                    std::move(single_result));
                std::tie(cont, turning, scan_value) = get_scan_val();
            }
            result_object[key] = std::move(result_unit);
            // reset initial guess
            omega = omega_initial_guess;
        }
    }

    timer.start_timing("Output");
    std::ofstream output(output_filename);
    output << result.dump();
    timer.pause_timing("Output");

    timer.pause_timing("All");
    std::cout << '\n';
    timer.print();
    std::cout << '\n';
#endif
    return 0;
}
