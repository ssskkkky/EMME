#ifndef TIMER_H  // Replace MATRIX_H with your unique guard macro name
#define TIMER_H

#include <chrono>
#include <string>
#include <unordered_map>
#include <utility>

class Timer {
   public:
    static Timer& get_timer();
    void start_timing(std::string func_name);

    void pause_timing(std::string func_name);

    void reset(std::string func_name);

    void print();

    std::unordered_map<
        std::string,
        std::pair<decltype(std::chrono::high_resolution_clock::now() -
                           std::chrono::high_resolution_clock::now()),
                  decltype(std::chrono::high_resolution_clock::now())>>
        time_consuming;

   private:
    Timer();
    Timer(const Timer&) = delete;
    Timer(Timer&&) = delete;
    decltype(auto) operator=(const Timer&) = delete;
    decltype(auto) operator=(Timer&&) = delete;
};

#endif
