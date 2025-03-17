#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "Matrix.h"

// Function to read contents of a file line by line
std::string readInputFromFile(const std::string& filename);

namespace util {

/**
 * @brief Casting integral type to double, keeping float-point type unchanged
 *
 * @tparam T
 */
template <typename T>
using numeric_t = std::conditional_t<std::is_integral_v<T>, double, T>;

template <bool B,
          template <typename...>
          class TrueTemplate,
          template <typename...>
          class FalseTemplate,
          typename... Args>
struct lazy_conditional;

template <template <typename...> class TrueTemplate,
          template <typename...>
          class FalseTemplate,
          typename... Args>
struct lazy_conditional<true, TrueTemplate, FalseTemplate, Args...> {
    using type = TrueTemplate<Args...>;
};

template <template <typename...> class TrueTemplate,
          template <typename...>
          class FalseTemplate,
          typename... Args>
struct lazy_conditional<false, TrueTemplate, FalseTemplate, Args...> {
    using type = FalseTemplate<Args...>;
};
template <typename T>
using get_value_type = typename T::value_type;

template <typename T>
using get_float_t = typename lazy_conditional<std::is_scalar_v<T>,
                                              numeric_t,
                                              get_value_type,
                                              T>::type;

// auxiliary function for Gauss-Kronrod quadrature
namespace detail {

/**
 * @brief The storage for abscissa and weight of order N Gauss-Kronrod
 * integration method. The data for N=5, 15 is precomputed.
 *
 */
template <size_t N>
struct gauss_kronrod_detail {
    constexpr static std::array<double, (N + 1) / 2> abscissa();
    constexpr static std::array<double, ((N - 1) / 2 + 1) / 2> gauss_weight();
    constexpr static std::array<double, (N + 1) / 2> kronrod_weight();
};

template <>
struct gauss_kronrod_detail<5> {
    constexpr static std::array<double, 3> abscissa() {
        return {0., 0.57735026918962576, 0.92582009977255146};
    }

    constexpr static std::array<double, 1> gauss_weight() { return {1.}; }

    constexpr static std::array<double, 3> kronrod_weight() {
        return {
            0.62222222222222222,
            0.49090909090909091,
            0.19797979797979798,
        };
    }
};

template <>
struct gauss_kronrod_detail<15> {
    constexpr static std::array<double, 8> abscissa() {
        return {
            0.,
            0.20778495500789847,
            0.40584515137739717,
            0.58608723546769113,
            0.74153118559939444,
            0.86486442335976907,
            0.94910791234275852,
            0.99145537112081264,
        };
    }
    constexpr static std::array<double, 4> gauss_weight() {
        return {
            0.41795918367346939,
            0.38183005050511894,
            0.27970539148927667,
            0.12948496616886969,

        };
    }
    constexpr static std::array<double, 8> kronrod_weight() {
        return {
            2.09482141084727828e-01, 2.04432940075298892e-01,
            1.90350578064785410e-01, 1.69004726639267903e-01,
            1.40653259715525919e-01, 1.04790010322250184e-01,
            6.30920926299785533e-02, 2.29353220105292250e-02,
        };
    }
};

template <>
struct gauss_kronrod_detail<31> {
    constexpr static std::array<double, 16> abscissa() {
        return {

            0.0,
            0.1011420669187175,
            0.20119409399743452,
            0.29918000715316881,
            0.39415134707756337,
            0.48508186364023968,
            0.57097217260853885,
            0.65099674129741697,
            0.72441773136017005,
            0.79041850144246593,
            0.84820658341042722,
            0.8972645323440819,
            0.9372733924007059,
            0.96773907567913913,
            0.98799251802048543,
            0.99800229869339706

        };
    }
    constexpr static std::array<double, 8> gauss_weight() {
        return {0.20257824192556112, 0.19843148532711152, 0.18616100001556193,
                0.1662692058169939,  0.1395706779261542,  0.10715922046717143,
                0.07036604748810768, 0.030753241996119};
    }
    constexpr static std::array<double, 16> kronrod_weight() {
        return {
            0.10133000701479155,   0.100769845523875595,  0.099173598721791959,
            0.0966427269836236785, 0.093126598170825321,  0.0885644430562117706,
            0.083080502823133021,  0.0768496807577203789, 0.069854121318728259,
            0.0620095678006706403, 0.053481524690928087,  0.0445897513247648766,
            0.035346360791375846,  0.0254608473267153202, 0.0150079473293161225,
            0.00537747987292334899};
    }
};

/**
 * @brief Order N Gauss-Kronrod quadrature, with embedded Gauss quadrature order
 * = (N-1)/2
 *
 */
template <size_t N, typename Tx>
struct gauss_kronrod : gauss_kronrod_detail<N> {
    using base = gauss_kronrod_detail<N>;

    /**
     * @brief Core function of gauss-kronrod integrate method, it integrates the
     * given function on [-1, 1].
     *
     * @tparam Func Function type
     * @param func Integrand
     * @return a pair of integral and err
     */
    template <typename Func>
    auto static gauss_kronrod_basic(const Func& func)
        -> std::pair<decltype(std::declval<Func>()(std::declval<Tx>())), Tx> {
        auto f0 = func(Tx{});
        using Ty = decltype(f0);
        using Tc = get_float_t<Ty>;
        constexpr auto gauss_order = (N - 1) / 2;

        Ty gauss_integral = gauss_order & 1
                                ? static_cast<Tc>(base::gauss_weight()[0]) * f0
                                : Ty{};
        Ty kronrod_integral = static_cast<Tc>(base::kronrod_weight()[0]) * f0;

        for (size_t i = 1; i < base::abscissa().size(); ++i) {
            Ty f = func(base::abscissa()[i]) + func(-base::abscissa()[i]);
            gauss_integral +=
                (gauss_order - i) & 1
                    ? static_cast<Tc>(base::gauss_weight()[i / 2]) * f
                    : Ty{};
            kronrod_integral += static_cast<Tc>(base::kronrod_weight()[i]) * f;
        }

        return std::make_pair(
            kronrod_integral,
            std::max(
                static_cast<Tx>(std::abs(kronrod_integral - gauss_integral)),
                static_cast<Tx>(std::abs(kronrod_integral) *
                                std::numeric_limits<Tx>::epsilon() * 2)));
    }

    template <typename Func>
    auto static gauss_kronrod_adaptive(const Func& func,
                                       Tx a,
                                       Tx b,
                                       size_t max_subdivide,
                                       Tx abs_tol,
                                       Tx global_rel_tol,
                                       Tx precision_goal) {
        // call stack
        std::vector<std::array<Tx, 2>> pending_intervals;
        // quadrature sum
        decltype(std::declval<Func>()(std::declval<Tx>())) sum{};

        Tx inv_scale = 2. / (b - a);
        pending_intervals.push_back({a, b});
        while (!pending_intervals.empty()) {
            auto [l, r] = pending_intervals.back();
            pending_intervals.pop_back();

            const Tx mid = (r + l) / 2;
            const Tx scale = (r - l) / 2;
            auto normalize_func = [&](Tx x) { return func(scale * x + mid); };
            auto result = gauss_kronrod_basic(normalize_func);
            auto integral = result.first * scale;
            auto err = result.second * scale;

            if (std::fpclassify(abs_tol) == FP_ZERO) {
                abs_tol = std::abs(global_rel_tol * integral);
            }
            if (std::ldexp(scale, max_subdivide) > 0.99 * (b - a) &&
                err > abs_tol * inv_scale + precision_goal &&
                err > std::abs(global_rel_tol * integral) + precision_goal) {
                pending_intervals.push_back({mid, r});
                pending_intervals.push_back({l, mid});
            } else {
                sum += integral;
            }
        }

        return sum;
    }
};

}  // namespace detail

template <typename Func, typename Ta, typename Tb, typename Te>
auto integrate(const Func& func,
               Ta a,
               Tb b,
               Te tol,
               Te prec,
               std::size_t max_subdivide,
               std::size_t integration_start_points) {
    using n_type = typename std::common_type<numeric_t<Ta>, numeric_t<Tb>,
                                             numeric_t<Te>>::type;
    if (integration_start_points == 15) {
        using impl = detail::gauss_kronrod<15, n_type>;
        if (std::fpclassify(a) == FP_INFINITE ||
            std::fpclassify(b) == FP_INFINITE) {
            // either of the endpoints is infinity
            return impl::gauss_kronrod_adaptive(
                [&](n_type x) {
                    const n_type c = std::cos(x);
                    return func(std::tan(x)) / (c * c);
                },
                std::atan(static_cast<n_type>(a)),
                std::atan(static_cast<n_type>(b)), max_subdivide, n_type{}, tol,
                prec);
        } else {
            return impl::gauss_kronrod_adaptive(func, a, b, max_subdivide,
                                                n_type{}, tol, prec);
        }
    } else if (integration_start_points == 31) {
        using impl = detail::gauss_kronrod<31, n_type>;
        if (std::fpclassify(a) == FP_INFINITE ||
            std::fpclassify(b) == FP_INFINITE) {
            // either of the endpoints is infinity
            return impl::gauss_kronrod_adaptive(
                [&](n_type x) {
                    const n_type c = std::cos(x);
                    return func(std::tan(x)) / (c * c);
                },
                std::atan(static_cast<n_type>(a)),
                std::atan(static_cast<n_type>(b)), max_subdivide, n_type{}, tol,
                prec);
        } else {
            return impl::gauss_kronrod_adaptive(func, a, b, max_subdivide,
                                                n_type{}, tol, prec);
        }
    }
    throw std::runtime_error("integration_start_points should be 15 or 31");
    return std::complex{0.0, 0.0};  // we will never reach here
}

template <typename Func, typename Te>
auto integrate(const Func& func,
               Te tol,
               Te prec,
               std::size_t max_subdivide,
               std::size_t integration_start_points) {
    using n_type = Te;
    if (integration_start_points == 15) {
        using impl = detail::gauss_kronrod<15, n_type>;
        return impl::gauss_kronrod_adaptive(
            [&](n_type x) {
                const n_type c = std::cos(x);
                return func(std::tan(x)) / (c * c);
            },
            0, std::numbers::pi / 2.0, max_subdivide, n_type{}, tol, prec);
    } else if (integration_start_points == 31) {
        using impl = detail::gauss_kronrod<31, n_type>;
        return impl::gauss_kronrod_adaptive(
            [&](n_type x) {
                const n_type c = std::cos(x);
                return func(std::tan(x)) / (c * c);
            },
            0, std::numbers::pi / 2.0, max_subdivide, n_type{}, tol, prec);
    }
    throw std::runtime_error("integration_start_points should be 15 or 31");
    return std::complex{0.0, 0.0};  // we will never reach here
}

template <typename Func, typename Tx>
auto integrate_coarse(const Func& func,
                      Tx a,
                      Tx b,
                      Tx tol = std::sqrt(std::numeric_limits<Tx>::epsilon()),
                      size_t max_subdivide = 3)
    -> decltype(std::declval<Func>()(std::declval<Tx>())) {
    using impl = detail::gauss_kronrod<5, Tx>;

    return impl::gauss_kronrod_adaptive(func, a, b, max_subdivide, Tx{}, tol);
}

// NOTE: When calculating log of bessel_i, the branch is not determined
template <typename T>
auto bessel_i_helper(const T& z, bool log = false) {
    constexpr double THRESHOLD = 2.e+7;
    int n = std::floor(std::abs(z)) + 1;
    T p0 = 0.;
    T p1 = 1.;
    T p_tmp{};
    double test_1 = std::max(
        std::sqrt(THRESHOLD * std::abs(p1) * std::abs(p0 - 2.0 * n / z * p1)),
        THRESHOLD);

    for (; std::abs(p1) <= test_1; ++n) {
        p_tmp = p0 - 2.0 * n / z * p1;
        p0 = p1;
        p1 = p_tmp;
    }

    T y0{1.0 / p1}, y1{}, y_tmp;
    T mu = 0.;
    for (n--; n > 0; --n) {
        y_tmp = 2. * n / z * y0 + y1;
        y1 = y0;
        y0 = y_tmp;
        mu += 2. * (std::real(z) < 0 ? 1 - 2 * (n & 1) : 1) * y1;
    }

    if (log) {
        return std::array<T, 2>{
            std::log(y0 / (mu + y0)) - (std::real(z) < 0 ? z : -z),
            std::log(y1 / (mu + y0)) - (std::real(z) < 0 ? z : -z)};
    }
    mu = std::exp(std::real(z) < 0 ? z : -z) * (mu + y0);
    return std::array<T, 2>{y0 / mu, y1 / mu};
}

template <typename T>
auto bessel_i_alter_helper(const T& z) {
    constexpr double THRESHOLD = 2.e+7;
    int n = std::floor(std::abs(z)) + 1;
    T p0 = 0.;
    T p1 = 1.;
    T p_tmp{};
    double test_1 = std::max(
        std::sqrt(THRESHOLD * std::abs(p1) * std::abs(p0 - 2.0 * n / z * p1)),
        THRESHOLD);

    for (; std::abs(p1) <= test_1; ++n) {
        p_tmp = p0 - 2.0 * n / z * p1;
        p0 = p1;
        p1 = p_tmp;
    }

    T y0{1.0 / p1}, y1{}, y_tmp;
    T mu = 0.;
    for (n--; n > 0; --n) {
        y_tmp = 2. * n / z * y0 + y1;
        y1 = y0;
        y0 = y_tmp;
        mu += 2. * (std::real(z) < 0 ? 1 - 2 * (n & 1) : 1) * y1;
    }

    return std::array<T, 4>{y0, y1, mu + y0, std::real(z) < 0 ? z : -z};
}

// NOTE: When calculating log of bessel_j, the branch is not determined
template <typename T>
auto bessel_j_helper(const T& z, bool log = false) {
    constexpr double THRESHOLD = 2.e+7;
    int n = std::floor(std::abs(z)) + 1;
    T p0 = 0.;
    T p1 = 1.;
    T p_tmp{};
    double test_1 = std::max(
        std::sqrt(THRESHOLD * std::abs(p1) * std::abs(p0 - 2.0 * n / z * p1)),
        THRESHOLD);

    for (; std::abs(p1) <= test_1; ++n) {
        p_tmp = 2.0 * n / z * p1 - p0;
        p0 = p1;
        p1 = p_tmp;
    }
    using namespace std::complex_literals;
    std::array<std::complex<double>, 4> ip{1., 1.i, -1., -1.i};

    T y0{1.0 / p1}, y1{}, y_tmp;
    T mu = 0.;
    for (n--; n > 0; --n) {
        y_tmp = 2. * n / z * y0 - y1;
        y1 = y0;
        y0 = y_tmp;
        mu += 2. * ip[std::imag(z) < 0 ? n & 3 : (4 - n & 3) & 3] * y1;
    }

    if (log) {
        mu = 1.i * (std::imag(z) < 0 ? -z : z) + std::log(mu + y0);
        return std::array<T, 2>{std::log(y0) - mu, std::log(y1) - mu};
    }
    mu = std::exp(1.i * (std::imag(z) < 0 ? -z : z)) * (mu + y0);
    return std::array<T, 2>{y0 / mu, y1 / mu};
}

std::string get_date_string();

template <typename T, std::size_t N>
auto unpack(const std::array<T, N>& arr) {
    return ([&]<std::size_t... idx>(std::index_sequence<idx...>) {
        return std::tuple<const decltype(arr[idx])...>(arr[idx]...);
    })(std::make_index_sequence<N>{});
}

}  // namespace util

template <typename T>
std::ostream& operator<<(std::ostream& outputStream,
                         const std::vector<std::complex<T>>& vectorObj) {
    // Print all the elements of vector using loop
    for (auto elem : vectorObj) {
        outputStream << elem.real() << " " << elem.imag() << " ";
    }
    return outputStream;
}

template <typename T, typename A>
std::ostream& operator<<(std::ostream& output_stream,
                         const Matrix<std::complex<T>, A>& matrix) {
    for (unsigned int i = 0; i < matrix.getRows(); i++) {
        for (unsigned int j = 0; j < matrix.getCols(); j++) {
            output_stream << matrix(i, j).real() << " " << matrix(i, j).imag()
                          << " ";  // Separate real and imaginary parts
        }
        output_stream << std::endl;  // New line after each row
    }

    return output_stream;
}

#endif
