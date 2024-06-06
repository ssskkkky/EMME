#include "Timer.h"

#include <iomanip>
#include <iostream>

using namespace std::chrono;

Timer::Timer() {
    time_consuming["integration"] = std::make_pair(
        high_resolution_clock::duration::zero(), high_resolution_clock::now());

    time_consuming["linear solver"] = std::make_pair(
        high_resolution_clock::duration::zero(), high_resolution_clock::now());
}

void Timer::start_timing(std::string func_name) {
    if (!time_consuming.contains(func_name)) { entries.push_back(func_name); }
    time_consuming[func_name].second = high_resolution_clock::now();
}

void Timer::pause_timing(std::string func_name) {
    auto end_time = high_resolution_clock::now();
    auto elapsed_time = end_time - time_consuming[func_name].second;
    time_consuming[func_name].first += elapsed_time;
}

void Timer::reset(std::string func_name) {
    time_consuming[func_name].first =
        decltype(time_consuming[func_name].first){};
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
