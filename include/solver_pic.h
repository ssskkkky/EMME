#ifndef SOLVER_PIC_H
#define SOLVER_PIC_H

#include <array>
#include <complex>
#include <iostream>
#include <vector>

template <typename T>
struct Marker {
    using value_type = T;

    value_type x;
    value_type v;
};

template <typename T>
struct PIC_State {
    using value_type = T;
    using velocity_type = std::vector<value_type>;
    using marker_container_type = std::vector<Marker<value_type>>;
    using field_type = std::vector<std::complex<value_type>>;

    PIC_State(std::size_t marker_num, std::size_t field_grid_num)
        : markers(initial_marker()), field(initial_field()) {}

    void put_velocity(velocity_type& vs) {}

    template <typename U>
    void update(U&& velocity, value_type dt) {}

   private:
    marker_container_type initial_marker() {}
    field_type initial_field() {}

    marker_container_type markers;
    field_type field;
};

template <typename T>
struct Integrator {
    using state_type = T;
    using value_type = typename state_type::value_type;
    using velocity_type = typename state_type::velocity_type;
    static constexpr std::size_t order = 3;

    Integrator(state_type& initial_state,
               value_type upper_err_bound = 1.e-7,
               value_type lower_err_bound = 1.e-10)
        : current_dt(0.1),
          state(initial_state),
          upper_err_bound_(upper_err_bound),
          lower_err_bound_(lower_err_bound) {}

    void step(value_type dt) {
        ([&]<auto... p_idx>(std::index_sequence<p_idx...>) {
            (([&]<auto... k_idx>(std::index_sequence<k_idx...>) {
                 constexpr auto p = p_idx;
                 state.put_velocity(intermediates[p]);
                 state.update((... + (coef[p][k_idx] * intermediates[k_idx])),
                              coef[p][p + 1] * dt);
             })(std::make_index_sequence<p_idx + 1>{}),
             ...);
        })(std::make_index_sequence<order>{});
    }

    auto step_adaptive() {
        auto state_current = state;

        value_type err{};
        while (true) {
            step(current_dt);
            err = ([&]<auto... k_idx>(std::index_sequence<k_idx...>) {
                return state.get_update_err(
                    (... + (coef[order][k_idx] * intermediates[k_idx])),
                    current_dt);
            })(std::make_index_sequence<3>{});
            if (err < upper_err_bound_) { break; }
            current_dt *= .5;
            state = state_current;
        }

        auto step_dt = current_dt;
        if (err < lower_err_bound_) { current_dt *= 2.; }

        return step_dt;
    }

   private:
    value_type current_dt;
    state_type& state;
    value_type upper_err_bound_;
    value_type lower_err_bound_;
    std::array<velocity_type, 3> intermediates;

    static inline constexpr std::array<std::array<double, 4>, 4> coef{
        {{1, 0.62653829327080},
         {0, 1, -0.55111240553326},
         {0, 1.5220585509963, -0.52205855099628, 0.92457411226246},
         {1., 0.13686116839369, -1.1368611683937}}};
};

#endif  // SOLVER_PIC_H
