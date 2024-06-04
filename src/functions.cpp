#include "functions.h"

#include <ctime>
#include <iomanip>

namespace util {

std::string get_date_string() {
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);

    std::stringstream ss;
    ss << std::put_time(&tm, "%FT%T%z");
    auto time_string = ss.str();
    char c = time_string[time_string.size() - 5];
    if (c == '+' || c == '-') {
        time_string.insert(time_string.size() - 2, ":");
    }
    return time_string;
}

}  // namespace util
