#include "Timer.h"

#include <iostream>

using namespace std::chrono;

Timer::Timer() {
    time_consuming["integration"] = std::make_pair(
        high_resolution_clock::duration::zero(), high_resolution_clock::now());

    time_consuming["linear solver"] = std::make_pair(
        high_resolution_clock::duration::zero(), high_resolution_clock::now());
}

void Timer::start_timing(std::string func_name) {
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
    for (auto& [timer, time_pair] : time_consuming) {
        std::cout << timer << ":"
                  << duration<double, seconds::period>(time_pair.first).count()
                  << "s" << std::endl;
    }
}

Timer& Timer::get_timer() {
    static Timer timer{};
    return timer;
}
