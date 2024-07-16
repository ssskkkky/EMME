#ifndef SOLVER_PIC_H
#define SOLVER_PIC_H

#include <array>
#include <complex>
#include <iostream>
#include <random>
#include <vector>

#include "Arithmetics.h"
#include "DedicatedThreadPool.h"
#include "Parameters.h"

template <typename T>
struct PIC_State {
    using value_type = T;
    using complex_type = std::complex<value_type>;

    struct Marker {
        value_type eta;
        const value_type v_para;
        const value_type v_perp;
        complex_type w;

        Marker& operator=(const Marker& other) {
            eta = other.eta;
            w = other.w;
            return *this;
        }
    };

    struct MarkerExtra {
        const value_type velocity_dependence_of_magnetic_drift_frequency;
        const value_type diamagnetic_drift_frequency;
        value_type bessel_j0;
        complex_type drift_center_pull_back;

        MarkerExtra& operator=(const MarkerExtra& other) {
            bessel_j0 = other.bessel_j0;
            drift_center_pull_back = other.drift_center_pull_back;
            return *this;
        }
    };

    using marker_container_type = std::vector<Marker>;
    using extra_container_type = std::vector<MarkerExtra>;
    using field_type = std::vector<complex_type>;

    struct velocity_type : util::ExpressionTemplate {
        std::vector<complex_type> data;

        velocity_type(std::size_t n) : data(n) {}

        decltype(auto) operator[](std::size_t idx) const { return data[idx]; }
        decltype(auto) operator[](std::size_t idx) { return data[idx]; }
    };

    PIC_State(const Parameters& para_input, std::size_t marker_num_per_cell)
        : para(para_input),
          cell_width(2 * para.length / para.npoints),
          markers(initialize_marker(marker_num_per_cell * para.npoints)),
          marker_extras(initialize_marker_extras()),
          qn_coef(cal_qn_coef()),
          field(initialize_field(para.npoints)) {}

    PIC_State& operator=(const PIC_State& other) {
        markers = other.markers;
        marker_extras = other.marker_extras;
        field = other.field;
        return *this;
    }

    inline auto initial_velocity_storage() const {
        return velocity_type(marker_num());
    }

    void put_velocity(velocity_type& vs) const {
        auto& thread_pool = DedicatedThreadPool<void>::get_instance();

        auto cal_velocity = [this, &vs](std::size_t begin, std::size_t end) {
            for (std::size_t i = begin; i < end; ++i) {
                const auto& [eta, v_para, v_perp, w] = markers[i];

                const auto x_perp = v_perp / para.vt;
                const auto sb = std::sqrt(para.b_theta *
                                          (1. + std::pow(para.shat * eta, 2)));
                const auto dj0 = -para.b_theta * para.shat * para.shat *
                                 x_perp * eta *
                                 std::cyl_bessel_j(1, x_perp * sb) / sb;

                const auto [cell_idx, cell_w] = locate(eta);
                const auto nf = field.size();
                const auto phi = (1. - cell_w) * field[cell_idx] +
                                 cell_w * field[(cell_idx + 1) % nf];
                const auto dphi =
                    ((1. - cell_w) * (field[(cell_idx + 1) % nf] -
                                      field[(cell_idx + nf - 1) % nf]) +
                     cell_w * (field[(cell_idx + 2) % nf] - field[cell_idx])) /
                    (2. * cell_width);

                const auto& [omega_dv, omega_st, j0, dc_pb] = marker_extras[i];
                vs[i] = std::conj(dc_pb) *
                        (complex_type{0., 1.} *
                             (omega_st - omega_d(eta) * omega_dv) * j0 * phi -
                         v_para / (para.q * para.R) * (j0 * dphi + dj0 * phi));
            }
        };

        constexpr std::size_t block_size = 512;
        std::vector<std::future<void>> res;
        for (std::size_t i = 0; i < marker_num() / block_size; ++i) {
            res.push_back(thread_pool.queue_task([&, i]() {
                cal_velocity(i * block_size, (i + 1) * block_size);
            }));
        }

        cal_velocity(marker_num() / block_size * block_size, marker_num());
        for (auto& f : res) { f.get(); }
    }

    template <typename U>
    void update(U&& velocity, value_type dt) {
        for (std::size_t i = 0; i < marker_num(); ++i) {
            auto& eta = markers[i].eta;
            eta = bound(eta + markers[i].v_para * dt / (para.q * para.R));
            markers[i].w += velocity[i] * dt;
        }

        solve_field();
    }

    template <typename U>
    auto get_update_err(U&& velocity, value_type dt) {
        value_type err{};
        value_type total{};
        for (std::size_t i = 0; i < field.size(); ++i) {
            err += std::real(velocity[i] * dt * std::conj(velocity[i] * dt));
            total += std::real(markers[i].w * std::conj(markers[i].w));
        }
        return std::sqrt(err / total);
    }

    // properties

    auto marker_num() const noexcept { return markers.size(); }

    // diagnostic

    decltype(auto) current_field() const { return field; }

   private:
    auto initialize_marker(std::size_t n) {
        marker_container_type initial_markers;
        initial_markers.reserve(n);

        auto&& gen = random_gen();

        std::uniform_real_distribution<value_type> uniform_eta(-para.length,
                                                               para.length);
        std::uniform_real_distribution<value_type> uniform_w(0, 0.001);
        std::normal_distribution<value_type> normal(0, para.vt);
        std::chi_squared_distribution<value_type> chi(2);

        for (std::size_t idx = 0; idx < n; ++idx) {
            initial_markers.push_back({.eta = uniform_eta(gen),
                                       .v_para = normal(gen),
                                       .v_perp = para.vt * std::sqrt(chi(gen)),
                                       .w = uniform_w(gen)});
        }

        return initial_markers;
    }
    auto initialize_marker_extras() {
        extra_container_type initial_extrasa;
        initial_extrasa.reserve(marker_num());
        for (const auto& [eta, v_para, v_perp, w] : markers) {
            initial_extrasa.push_back(
                {.velocity_dependence_of_magnetic_drift_frequency =
                     (v_para * v_para + .5 * v_perp * v_perp) /
                     (2. * para.vt * para.vt),
                 .diamagnetic_drift_frequency =
                     para.omega_s_i *
                     (1. + para.eta_i * ((v_para * v_para + v_perp * v_perp) /
                                             (2. * para.vt * para.vt) -
                                         1.5))});
        }
        return initial_extrasa;
    }
    auto initialize_field(std::size_t n) {
        field_type initial_field(n);

        return initial_field;
    }

    auto inline locate(auto eta) const {
        std::size_t idx = (eta + para.length) / cell_width;
        value_type w = (eta + para.length) / cell_width - idx;
        return std::make_pair(idx, w);
    }

    void solve_field() {
        constexpr std::size_t batch_count = 1 << 8;
        auto& thread_pool = DedicatedThreadPool<void>::get_instance();

        const auto nf = field.size();
        // buffer
        static std::vector<complex_type> buffer(batch_count * nf);

        // calculate density
        auto cal_density = [this](std::size_t begin, std::size_t end,
                                  std::size_t buffer_begin) {
            for (std::size_t i = 0; i < field.size(); ++i) {
                buffer[buffer_begin + i] = 0;
            }
            for (std::size_t i = begin; i < end; ++i) {
                const auto& [eta, v_para, v_perp, w] = markers[i];
                auto& [omega_dv, omega_st, j0, dc_pb] = marker_extras[i];

                const auto x_perp = v_perp / para.vt;
                const auto sb = std::sqrt(para.b_theta *
                                          (1. + std::pow(para.shat * eta, 2)));
                j0 = std::cyl_bessel_j(0, x_perp * sb);
                dc_pb = std::exp(complex_type{
                    0., -omega_d_integral(eta, v_para) * omega_dv});

                const auto den = j0 * w * dc_pb;

                const auto [cell_idx, cell_w] = locate(eta);

                // left grid point
                buffer[buffer_begin + cell_idx] += den * (1. - cell_w);
                buffer[buffer_begin + (cell_idx + 1) % field.size()] +=
                    den * cell_w;
            }
        };

        std::vector<std::future<void>> res;
        const auto batch_size = marker_num() / batch_count;
        const auto remain_size = marker_num() % batch_count;
        for (std::size_t i = 0; i < batch_count; ++i) {
            res.push_back(thread_pool.queue_task([&, i]() {
                cal_density(
                    i * batch_size + (i < remain_size ? i : remain_size),
                    (i + 1) * batch_size +
                        (i < remain_size ? i + 1 : remain_size),
                    i * nf);
            }));
        }
#if 0
        constexpr std::size_t reduce_size =
            4;  // batch_count should be a power of this number
        // std::vector<std::future<void>> reduce_res;
        // NOTE: Reservation is necessary here since I do not lock this vector
        // when visiting it in other threads
        res.reserve((batch_count * reduce_size - 1) / (reduce_size - 1));
        for (std::size_t i = 0, current_level_idx = 0,
                         current_level_task_count = batch_count / reduce_size;
             i < (batch_count - 1) / (reduce_size - 1);
             ++i, ++current_level_idx) {
            if (current_level_idx == current_level_task_count) {
                current_level_task_count /= reduce_size;
                current_level_idx = 0;
            }
            res.push_back(thread_pool.queue_task([&, current_level_idx,
                                                  current_level_task_count]() {
                // wait for previous tasks
                std::size_t last_level_begin =
                    (batch_count - current_level_task_count * reduce_size) *
                    reduce_size / (reduce_size - 1);
                for (std::size_t j = 0; j < reduce_size; ++j) {
                    res[last_level_begin + current_level_idx * reduce_size + j]
                        .get();
                }
                std::size_t block_length =
                    batch_count / current_level_task_count * nf;
                for (std::size_t j = 1; j < reduce_size; ++j) {
                    for (std::size_t k = 0; k < nf; ++k) {
                        buffer[current_level_idx * block_length + k] +=
                            buffer[current_level_idx * block_length +
                                   j * block_length / reduce_size + k];
                    }
                }
            }));
        }

        res.back().get();
        for (std::size_t idx = 0; idx < nf; ++idx) {
            field[idx] = buffer[idx] * qn_coef[idx];
        }
#else
        for (auto& f : res) { f.get(); }
        for (std::size_t idx = 0; idx < nf; ++idx) { field[idx] = 0; }
        for (std::size_t i = 0; i < batch_count; ++i) {
            for (std::size_t idx = 0; idx < nf; ++idx) {
                field[idx] += buffer[i * nf + idx];
            }
        }

        for (std::size_t idx = 0; idx < nf; ++idx) {
            field[idx] *= qn_coef[idx];
        }
#endif
    }

    static decltype(auto) random_gen() {
        static std::mt19937 gen{std::random_device{}()};
        return gen;
    }

    inline auto omega_d(auto eta) const {
        return 2. * para.epsilon_n * para.omega_s_i *
               (std::cos(eta) + para.shat * eta * std::sin(eta));
    }

    inline auto omega_d_integral(auto eta, auto v_para) const {
        return (para.q * para.R / v_para) * 2 * para.epsilon_n *
               para.omega_s_i *
               (std::sin(eta) * (1. + para.shat) -
                para.shat * eta * std::cos(eta));
    }

    inline auto cal_qn_coef() const {
        std::vector<value_type> gamma0(para.npoints);
        for (std::size_t idx = 0; idx < qn_coef.size(); ++idx) {
            const auto b =
                para.b_theta *
                (1. +
                 std::pow(para.shat * (idx * cell_width - para.length), 2));
            gamma0[idx] = std::cyl_bessel_i(0, b) * std::exp(-b);
            // (1+\frac{1}{\tau}+\Gamma_0)\delta\phi = \delta n
            // ---------------------------
            //             ||
            //             ||
            //          gamma0 now equals the inverse of this coefficient
            //          (with
            //              proper normalization)
            gamma0[idx] = 2. * para.length /
                          (marker_num() * (1. + 1. / para.tau - gamma0[idx]) *
                           cell_width);
        }
        return gamma0;
    }

    inline auto bound(auto eta) {
        eta = std::fmod(eta + para.length, 2 * para.length);
        return eta < 0 ? eta + para.length : eta - para.length;
    }

    const Parameters& para;
    const value_type cell_width;
    marker_container_type markers;
    extra_container_type marker_extras;
    const std::vector<value_type> qn_coef;
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
          lower_err_bound_(lower_err_bound),
          intermediates(([&]<auto... p_idx>(std::index_sequence<p_idx...>) {
              return std::array<velocity_type, order>{
                  ((void)p_idx, state.initial_velocity_storage())...};
          })(std::make_index_sequence<order>{})) {}

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
