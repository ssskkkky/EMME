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

/**
 * @brief 计算给定a值的函数 cos(x) + a*x*sin(x)
 * @param x 函数的输入变量x
 * @param a 函数的参数a (非负实数)
 * @return 函数在x处的值
 */
double omegad_potential_f(double x, double a) {
    return std::cos(x) + a * x * std::sin(x);
}

/**
 * @brief 使用二分法在 [0, pi] 区间内找到函数 cos(x) + a*x*sin(x) 的零点。
 * 函数假定在给定区间内有且只有一个零点。
 * @param a 函数的参数a (非负实数)
 * @param tolerance 寻找零点的精度容差
 * @param maxIterations 最大迭代次数，防止无限循环
 * @return 找到的零点近似值
 */
double findZeroPoint(double a, double tolerance, int maxIterations) {
    double low = 0.0;
    double high = M_PI;

    // 检查区间端点是否已经是零点（或非常接近）
    if (std::fabs(omegad_potential_f(low, a)) < tolerance) { return low; }
    if (std::fabs(omegad_potential_f(high, a)) < tolerance) { return high; }

    double mid = 0.0;
    for (int i = 0; i < maxIterations; ++i) {
        mid = low + (high - low) / 2.0;  // 更稳定的中点计算方式
        double f_mid = omegad_potential_f(mid, a);

        // 如果中点的值已经足够接近零，或者区间足够小，则认为找到了零点
        if (std::fabs(f_mid) < tolerance || (high - low) / 2.0 < tolerance) {
            return mid;
        }

        // 判断零点在哪一半区间
        // f(low) 和 f(high) 总是异号的 (f(0)=1, f(pi)=-1)
        // 如果 f(low) * f(mid) < 0，说明零点在 [low, mid] 区间
        if (omegad_potential_f(low, a) * f_mid < 0) {
            high = mid;
        } else {  // 否则，零点在 [mid, high] 区间
            low = mid;
        }
    }

    // 达到最大迭代次数后，返回当前最佳估计值
    return mid;
}

/**
 * @brief 计算函数 cos(x) + a*x*sin(x) 从 0 到其第一个零点的平均值。
 * @param a 函数的参数a (非负实数)
 * @return 函数在指定区间内的平均值
 */
double calculateAverageValue(double a) {
    double x0 = findZeroPoint(a);  // 首先找到第一个零点

    // 根据公式计算定积分的值： (1+a)*sin(x0) - a*x0*cos(x0)
    // F(0) 值为 0，所以只需要计算 F(x0)
    double integral_value = (1.0 + a) * std::sin(x0) - a * x0 * std::cos(x0);

    // 平均值 = 积分值 / (x0 - 0)
    // 注意：由于 f(0)=1 且 f(pi)=-1，且 a >= 0，零点 x0 总是 >
    // 0，所以不会出现除以零的情况。
    return integral_value / x0;
}
