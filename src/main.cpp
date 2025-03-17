#include <cfenv>
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

auto solve_once_eigen(const auto& input,
                      auto& omega_initial_guess,
                      std::ofstream& eigen_matrix_file) {
    auto& timer = Timer::get_timer();
    double tol = input.at("iteration_precision");

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

auto solve_once_pic(const auto& input,
                    auto&,
                    std::ofstream& eigen_matrix_file) {
    auto& timer = Timer::get_timer();
    timer.start_timing("Initial");

    auto& para = Parameters::generate(input);
    const std::size_t marker_per_cell = input.at("marker_per_cell");
    PIC_State<double> state(para, marker_per_cell);
    Integrator integrator(state);

    const std::size_t nt = input.at("step_number");
    const double dt = input.at("time_step");

    std::vector<std::array<double, 3>> stats;
    stats.reserve(nt);

    timer.pause_timing("Initial");
    for (std::size_t idx = 0; idx < nt; ++idx) {
        integrator.step(dt);

        // diagnostics
        timer.start_timing("Diagnostics");
        const auto& current_field = state.current_field();
        const auto nf = current_field.size();
        eigen_matrix_file.write(
            reinterpret_cast<const char*>(current_field.data()),
            sizeof(current_field[0]) * nf);

        auto [real, imag, norm] = std::accumulate(
            current_field.begin(), current_field.end(), std::array<double, 3>{},
            [](auto acc, const auto& val) {
                return std::array{acc[0] + std::real(val),
                                  acc[1] + std::imag(val),
                                  acc[2] + std::real(val * std::conj(val))};
            });
        stats.push_back({real / nf, imag / nf, std::sqrt(norm / nf)});
        std::cout << "        " << idx + 1 << '/' << nt
                  << " phi[0]: " << current_field[nf / 2] << '\n';
        timer.pause_timing("Diagnostics");
    }

    auto eigen_value = util::calculate_omega(stats, dt);
    std::cout << "        Eigenvalue: " << eigen_value << '\n';

    auto single_result = Value::create_object();
    auto& eva = single_result["eigenvalue"] = Value::create_array(2);

    eva[0] = eigen_value.real();
    eva[1] = eigen_value.imag();

    single_result["eigenvector"] =
        Value::create_typed_array(state.current_field());

    return single_result;
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
    auto input_all = util::json::parse_file(filename);

    auto invoke_solver = [&]<typename... Args>(Args&&... args) {
        std::string method = input_all.at("method");
        if ("eigen" == method) {
            return solve_once_eigen(std::forward<Args>(args)...);
        } else if ("PIC" == method) {
            return solve_once_pic(std::forward<Args>(args)...);
        }
        std::ostringstream oss;
        oss << "Method '" << input_all.at("method").as_string()
            << "' is not supported, yet.\n";
        throw std::runtime_error(oss.str());
    };

    auto& timer = Timer::get_timer();
    timer.start_timing("All");

    std::string output_filename = "output.json";

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

    // scan_config["key"] = {head, step, tail, another_tail};
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
            invoke_solver(input_all, omega_initial_guess, eigen_matrix_file));
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
                    if (new_omega.is_string()) {
                        // NaN
                        omega = omega_initial_guess;
                    } else {
                        omega.real(new_omega[0]);
                        omega.imag(new_omega[1]);
                    }
                }

                std::cout << "    " << key << ":" << scan_value << '\n';

                auto eigen_matrix_file_name = "eigenMatrics/" + key + "Eq" +
                                              std::to_string(scan_value) +
                                              ".bin";
                std::ofstream eigen_matrix_file(eigen_matrix_file_name,
                                                std::ios::binary);
                try {
                    auto single_result =
                        invoke_solver(input, omega, eigen_matrix_file);
                    single_result["eigenMatrix"] =
                        eigen_matrix_file
                            ? eigen_matrix_file_name
                            : "Can not open '" + eigen_matrix_file_name +
                                  "' for write.";
                    single_result["scan_value"] = scan_value;
                    scan_result_array.as_array().push_back(
                        std::move(single_result));
                } catch (const std::exception& e) {
                    auto err_result = Value::create_object();
                    err_result["eigenvalue"] = "NaN";
                    err_result["reason"] = e.what();
                    scan_result_array.as_array().push_back(
                        std::move(err_result));
                    std::cerr << "        " << e.what() << '\n';
                }
                std::tie(cont, turning, scan_value) = get_scan_val();
            }
            result_object[key] = std::move(result_unit);
            // reset initial guess
            omega = omega_initial_guess;
        }
    }

    timer.start_timing("Output");
    std::ofstream output(output_filename);
    output << result.pretty_print();
    timer.pause_timing("Output");

    timer.pause_timing("All");
    std::cout << '\n';
    timer.print();
    std::cout << '\n';

    return 0;
}
