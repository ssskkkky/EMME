#include <cmath>
#include <iostream>
#include <limits>

#include "solver_pic.h"

struct double_state {
    using value_type = double;

    // velocity type need to support value_type*velocity_type and
    // velocity_type + velocity_type
    // expression template can also fit in here, as long as `update` method
    // accepts it
    struct velocity_type {
        value_type v0, v1;
    };

    friend auto operator*(value_type a, velocity_type v) {
        return velocity_type{a * v.v0, a * v.v1};
    }
    friend auto operator+(velocity_type lhs, velocity_type rhs) {
        return velocity_type{lhs.v0 + rhs.v0, lhs.v1 + rhs.v1};
    }

    // d^2x/dt^2 = -x
    void put_velocity(velocity_type& v) {
        v.v0 = x1;
        v.v1 = -x0;
    }

    void update(velocity_type v, value_type dt) {
        x0 += v.v0 * dt;
        x1 += v.v1 * dt;
        t += dt;
    }

    auto get_update_err(velocity_type v, value_type dt) {
        auto l2 = [](value_type a, value_type b) {
            return std::sqrt(.5 * (a * a + b * b));
        };
        return l2(x0, x1) < std::numeric_limits<value_type>::epsilon()
                   ? l2(v.v0 * dt, v.v1 * dt)
                   : l2(v.v0 * dt, v.v1 * dt) / l2(x0, x1);
    }

    double t;
    double x0, x1;
};

int main() {
    constexpr double total_t = 10;
    {
        double_state x2{0, 0, 1};
        Integrator<double_state> f2(x2, 1.e-5, 1.e-7);

        constexpr double dt = 0.01;
        for (std::size_t i = 0; i < total_t / dt; ++i) {
            f2.step(dt);
            if (x2.x0 - std::sin(x2.t) > 1.e-5) {
                std::cout << x2.t << ", " << dt << ", ";
                std::cout << x2.x0 << "-> " << x2.x0 - std::sin(x2.t) << '\n';
            }
        }
    }
    {
        double_state x2{0, 0, 1};
        Integrator<double_state> f2(x2, 1.e-5, 1.e-7);

        std::size_t c = 0;
        while (x2.t < total_t) {
            auto dt = f2.step_adaptive();
            if (x2.x0 - std::sin(x2.t) > 1.e-5) {
                std::cout << x2.t << ", " << dt << ", ";
                std::cout << x2.x0 << "-> " << x2.x0 - std::sin(x2.t) << '\n';
            }
            ++c;
        }
        std::cout << "Adpative method use " << c << " steps.\n";
    }

    return 0;
}
