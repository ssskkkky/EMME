#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

class Timer {
   public:
    static Timer& get_timer();
    void start_timing(std::string func_name);
    void pause_and_start(std::string func_name);
    void pause_timing();
    void pause_timing(std::string func_name);
    void reset();
    void print();

   private:
    Timer() = default;
    Timer(const Timer&) = delete;
    Timer(Timer&&) = delete;
    decltype(auto) operator=(const Timer&) = delete;
    decltype(auto) operator=(Timer&&) = delete;

    std::vector<std::string> entries;
    std::unordered_map<
        std::string,
        std::pair<std::chrono::high_resolution_clock::duration,
                  std::chrono::high_resolution_clock::time_point>>
        time_consuming;

    std::string current_func_name;
};

#endif
