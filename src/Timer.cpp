#include "Timer.h"

#include <iomanip>
#include <iostream>

using namespace std::chrono;

void Timer::start_timing(std::string func_name) {
    current_func_name = func_name;
    auto emplace_reuslt = time_consuming.emplace(
        func_name, std::make_pair(high_resolution_clock::duration::zero(),
                                  high_resolution_clock::now()));
    if (emplace_reuslt.second) {
        entries.push_back(func_name);
    } else {
        emplace_reuslt.first->second.second = high_resolution_clock::now();
    }
}

void Timer::pause_and_start(std::string func_name) {
    pause_timing();
    start_timing(func_name);
}

void Timer::pause_timing() {
    auto end_time = high_resolution_clock::now();
    auto elapsed_time = end_time - time_consuming.at(current_func_name).second;
    time_consuming.at(current_func_name).first += elapsed_time;
}

void Timer::pause_timing(std::string func_name) {
    auto end_time = high_resolution_clock::now();
    auto elapsed_time = end_time - time_consuming.at(func_name).second;
    time_consuming.at(func_name).first += elapsed_time;
}

void Timer::reset() {
    entries.clear();
    time_consuming.clear();
}

void Timer::print() {
    std::size_t max_length = 0;
    for (auto& name : entries) {
        max_length = max_length < name.size() ? name.size() : max_length;
    }

    std::cout << '+';
    for (std::size_t i = 0; i < max_length + 16; ++i) { std::cout << '-'; }
    std::cout << "+\n"
              << std::left << std::setw(max_length + 17)
              << "| Time consumption";
    std::cout << "|\n+";
    for (std::size_t i = 0; i < max_length; ++i) { std::cout << "-"; }
    std::cout << "--+-------------+\n";
    for (const auto& name : entries) {
        auto& time_used = time_consuming.at(name).first;
        std::cout << std::left << "| " << std::setw(max_length) << name << " | "
                  << std::setw(11)
                  << duration<double, seconds::period>(time_used).count()
                  << "s|\n";
    }
    std::cout << '+';
    for (std::size_t i = 0; i < max_length; ++i) { std::cout << "-"; }
    std::cout << "--+-------------+\n";
}

Timer& Timer::get_timer() {
    static Timer timer{};
    return timer;
}
