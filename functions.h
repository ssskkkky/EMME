#ifndef FUNCTIONS_H  // Replace MATRIX_H with your unique guard macro name
#define FUNCTIONS_H

#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>

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
                                       Tx local_abs_tol,
                                       Tx global_rel_tol)
        -> decltype(std::declval<Func>()(std::declval<Tx>())) {
        const Tx mid = (b + a) / 2;
        const Tx scale = (b - a) / 2;
        auto normalize_func = [&](Tx x) { return func(scale * x + mid); };

        auto result = gauss_kronrod_basic(normalize_func);
        auto integral = result.first * scale;
        auto err = result.second * scale;
        if (std::fpclassify(local_abs_tol) == FP_ZERO) {
            local_abs_tol = std::abs(global_rel_tol * integral);
        }

        if (max_subdivide > 0 && err > local_abs_tol &&
            err > std::abs(global_rel_tol * integral)) {
            return gauss_kronrod_adaptive(func, a, mid, max_subdivide - 1,
                                          local_abs_tol / 2, global_rel_tol) +
                   gauss_kronrod_adaptive(func, mid, b, max_subdivide - 1,
                                          local_abs_tol / 2, global_rel_tol);
        }

#ifdef _TRACE
        std::cout << "[TRACE] Gauss-Kronrod subdivide stopped at [" << a << ", "
                  << b << "] with remaining subdivision = " << max_subdivide
                  << ", err = " << err << '\n';
#endif

        return integral;
    }
};

}  // namespace detail

template <typename Func, typename Ta, typename Tb, typename Te>
auto integrate(const Func& func,
               Ta a,
               Tb b,
               Te tol,
               size_t max_subdivide = 15) {
    using n_type = typename std::common_type<numeric_t<Ta>, numeric_t<Tb>,
                                             numeric_t<Te>>::type;
    using impl = detail::gauss_kronrod<15, n_type>;
    if (b + 1 == b) {
        // b is infinity
        return impl::gauss_kronrod_adaptive(
            [&](n_type x) {
                const n_type c = std::cos(x);
                return func(std::tan(x)) / (c * c);
            },
            std::atan(static_cast<n_type>(a)),
            std::atan(static_cast<n_type>(b)), max_subdivide, n_type{}, tol);
    } else {
        return impl::gauss_kronrod_adaptive(func, a, b, max_subdivide, n_type{},
                                            tol);
    }
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

/**
 * @brief Modified Bessel function of the first kind with integer order and
 * complex argument. The evaluation is done using contour integral.
 *
 * @param i order
 * @param z argument
 * @param eps accuracy, default to be 1e-7
 */
template <typename T>
auto bessel_ic(int i, T z, double eps = 1.e-7) {
    using f_type = get_float_t<T>;
    constexpr auto pi = static_cast<f_type>(M_PI);
    return integrate(
               [=](f_type theta) {
                   return std::exp(z * std::cos(theta)) * std::cos(i * theta);
               },
               0, pi, static_cast<f_type>(eps)) /
           pi;
}

}  // namespace util

#endif
