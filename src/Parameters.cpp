#include "Parameters.h"

#include <cmath>
#include <tuple>

#include "functions.h"
using namespace std::literals;

const Parameters& Parameters::generate(const util::json::Value& input) {
    // Stellarator is the largest derived class
    alignas(Stellarator) static std::byte buffer[sizeof(Stellarator)];
    auto para_ptr = reinterpret_cast<Parameters*>(buffer);

    // Parameters and Stellarator are both trivially
    // destructible, no need to bother calling their
    // destructors.
    if (std::string{"tokamak"}.compare(input.at("conf")) == 0) {
        new (buffer) Parameters(input);
    } else if (std::string{"stellarator"}.compare(input.at("conf")) == 0) {
        new (buffer) Stellarator(input);
    } else if (std::string{"cylinder"}.compare(input.at("conf")) == 0) {
        new (buffer) Cylinder(input);
    } else {
        throw std::runtime_error("Input configuration not supported yet.");
    }

    return *para_ptr;
}

Parameters::Parameters(const util::json::Value& input)
    : q(input.at("q")),
      shat(input.at("shat")),
      tau(input.at("tau")),
      epsilon_n(input.at("epsilon_n")),
      epsilon_r(input.at("epsilon_r")),
      eta_i(input.at("eta_i")),
      eta_e(input.at("eta_e")),
      b_theta(input.at("k_rho") * input.at("k_rho")),
      beta_e(input.at("beta_e")),
      R(input.at("R")),
      vt(input.at("vt")),
      omega_d_coeff(input.at("omega_d_coeff")),
      length(input.at("length")),
      theta(input.at("theta")),
      npoints(input.at("npoints")),
      iteration_step_limit(input.at("iteration_step_limit")),
      integration_precision(input.at("integration_precision")),
      integration_accuracy(input.at("integration_accuracy")),
      integration_iteration_limit(input.at("integration_iteration_limit")),
      integration_start_points(input.at("integration_start_points")),
      arc_coeff(input.at("arc_coeff")),
      alpha(q * q * R * beta_e / (epsilon_n * R) *
            ((1 + eta_e) + 1 / tau * (1 + eta_i))),
      water_bag_weight_vpara(input.at("water_bag_weight_vpara")),
      water_bag_weight_vperp(input.at("water_bag_weight_vperp")),
      omega_s_i(-(std::sqrt(b_theta) * vt) / (epsilon_n * R)),
      omega_s_e(-tau * omega_s_i),
      omega_d_bar(2.0 * epsilon_n * omega_s_i * omega_d_coeff),
      drift_center_transformation_switch(
          input.at("drift_center_transformation_switch").as_boolean()) {}

void Parameters::parameterInit() {
    alpha = q * q * R * beta_e / (epsilon_n * R) *
            ((1 + eta_e) + 1 / tau * (1 + eta_i));
    omega_s_i = -(std::sqrt(b_theta) * vt) / (epsilon_n * R);
    omega_s_e = -tau * omega_s_i;
    omega_d_bar = 2.0 * epsilon_n * omega_s_i * omega_d_coeff;
};

double Parameters::g_integration_f(double eta) const {
    return -((alpha * eta) / 2.0) + shat * theta * std::cos(eta) -
           shat * eta * std::cos(eta) + std::sin(eta) + shat * std::sin(eta) +
           0.25 * alpha * std::sin(2.0 * eta) -
           (1 - shat) * epsilon_r * eta / q / q;
    // One should check the magnetic effect in the last term
    // (1 - shat) * epsilon_r*eta / q / q;. Maybe add alpha somehow in this term
}

double Parameters::beta_1(double eta, double eta_p) const {
    return (q * R) / vt * (omega_d_bar) *
           (g_integration_f(eta) - g_integration_f(eta_p));
}

double Parameters::beta_1_e(double eta, double eta_p) const {
    return (q * R) / vt * (omega_d_bar * omega_s_e / omega_s_i) *
           (g_integration_f(eta) - g_integration_f(eta_p));
}

double Parameters::bi(double eta) const {
    return b_theta *
           (1.0 + pow(shat * (eta - theta) - alpha * std::sin(eta), 2));
}
std::complex<double> Parameters::lambda_f_tau(double eta,
                                              double eta_p,
                                              std::complex<double> tau) const {
    return 1.0 + 0.5 * std::complex<double>(0.0, 1.0) * (tau * vt) /
                     (q * R * (eta - eta_p)) * beta_1(eta, eta_p);
}

std::complex<double> Parameters::h_f_tau(std::complex<double> omega,
                                         std::complex<double> tau) const {
    return exp(std::complex<double>(0.0, 1.0) * tau * omega);
}

std::complex<double> Parameters::kappa_f_tau(unsigned int m,
                                             double eta,
                                             double eta_p,
                                             std::complex<double> omega) const

{
    // Define the integrand function
    auto integrand = [&](double taut_transformed) {
        const auto omi = -std::copysign(1, omega.real());
        const auto exp_arg =
            std::exp(-omi * 1.i * std::atan(taut_transformed / arc_coeff));
        const auto taut = taut_transformed * exp_arg;

        const auto jacob =
            exp_arg - (1.i * exp_arg * omi * taut_transformed) /
                          (arc_coeff *
                           (1.0 + std::pow((taut_transformed / arc_coeff), 2)));

        std::complex<double> lambda_f_tau_term = lambda_f_tau(eta, eta_p, taut);
        const auto bi_eta = bi(eta);
        const auto bi_eta_p = bi(eta_p);

        const auto [y0, y1, mu, z] = util::bessel_i_alter_helper(
            std::sqrt(bi_eta * bi_eta_p) / lambda_f_tau_term);

        const auto lambda_f_tau_term_cubic_inv =
            std::pow(lambda_f_tau_term, -3.);
        const auto norm_vel = (q * R * (eta - eta_p)) / (vt * taut);

        const std::complex<double> i0_coef =
            (omega -
             omega_s_i * (1.0 + eta_i * (0.5 * norm_vel * norm_vel - 1.5))) /
                lambda_f_tau_term +
            omega_s_i * eta_i * (.5 * (bi_eta + bi_eta_p) - lambda_f_tau_term) *
                lambda_f_tau_term_cubic_inv;

        const std::complex<double> i1_coef = -omega_s_i * eta_i *
                                             std::sqrt(bi_eta * bi_eta_p) *
                                             lambda_f_tau_term_cubic_inv;

        const auto beta_1_val = beta_1(eta, eta_p);

        // logarithmic of normalized parallel velocity in exponential term,
        // following terms is similar
        const auto log_norm_vel = -0.5 * norm_vel * norm_vel;
        const auto log_i_beta = -.5i * beta_1_val * norm_vel;
        const auto log_hf_tau = 1.i * taut * omega;
        const auto log_exp_term_int_lambda_tau =
            -(bi_eta + bi_eta_p) / (2.0 + 1.i * beta_1_val / norm_vel);

        const auto log_coef = log_norm_vel + log_i_beta + log_hf_tau +
                              log_exp_term_int_lambda_tau;

        // deal with underflow problem
        const auto safe_exp = [](auto var) {
            if (std::real(var) < -40.) {
                return std::complex<double>{0.};
            } else {
                return std::exp(var);
            }
        };
        return std::pow(norm_vel, m) / taut * jacob * safe_exp(log_coef - z) *
               (i0_coef * y0 + i1_coef * y1) / mu;
    };

    auto result =
        util::integrate(integrand, integration_precision, integration_accuracy,
                        integration_iteration_limit, integration_start_points);

    return -std::complex<double>(0, 1.0) * (q * R) /
           (vt * std::sqrt((2.0 * M_PI))) * result;
}

std::complex<double> Parameters::kappa_f_tau_e(unsigned int m,
                                               double eta,
                                               double eta_p,
                                               std::complex<double> omega) const

{
    switch (m) {
        case 0:
            return 0.0;
        case 1:
            return -std::complex<double>(0.0, 1.0) * (q * R) /
                   (2.0 * vt * tau) * (omega - omega_s_e) * (eta - eta_p) /
                   (std::abs(eta - eta_p));
        case 2:
            return (q * q * R * R) / (2.0 * vt * vt * tau) * (eta - eta_p) /
                   (std::abs(eta - eta_p)) *
                   (omega * (omega - omega_s_e) * (eta - eta_p) -
                    beta_1_e(eta, eta_p) * vt / (q * R) *
                        (omega - omega_s_e * (1.0 + eta_e)));
        default:
            // Handle unexpected mode (throw exception or return special value)
            throw std::invalid_argument("Unsupported mode value");
    }
}

Stellarator::Stellarator(const util::json::Value& input)
    : Parameters(input),
      eta_k(input.at("eta_k")),
      lh(input.at("lh")),
      mh(input.at("mh")),
      epsilon_h_t(input.at("epsilon_h_t")),
      alpha_0(input.at("alpha_0")),
      r_over_R(input.at("r_over_R")),
      deltap(-0.25 * alpha),
      beta_e_p(beta_e * (1.0 + eta_e) / (epsilon_n * R)),
      rdeltapp((-alpha + (2.0 * shat - 3) * deltap)),
      curvature_aver(mh / lh * r_over_R / (q * R) * (4.0 - shat) +
                     (-alpha + 2 * shat * deltap + 0) / R) {}

double Stellarator::sigma_f(double eta) const {
    return shat * (eta - eta_k) +
           (deltap * (1 + shat) + rdeltapp) * std::sin(eta);
}

double Stellarator::bi(double eta) const {
    return b_theta * (1.0 + pow(sigma_f(eta), 2));
}

void Stellarator::parameterInit() {
    alpha = q * q * R * beta_e / (epsilon_n * R) *
            ((1 + eta_e) + 1 / tau * (1 + eta_i));
    omega_s_i = -(std::sqrt(b_theta) * vt) / (epsilon_n * R);
    omega_s_e = -tau * omega_s_i;
    omega_d_bar = 2.0 * epsilon_n * omega_s_i * omega_d_coeff;

    deltap = -0.25 * alpha;
    beta_e_p = beta_e * (1.0 + eta_e) / (epsilon_n * R);
    rdeltapp = (-alpha + (2.0 * shat - 3) * deltap);
    curvature_aver = mh / lh * r_over_R / (q * R) * (4.0 - shat) +
                     (-alpha + 2 * shat * deltap + 0) / R;
};

double Stellarator::g_integration_f(double eta) const {
    return (eta * (-1 + lh - mh * q) * std::pow(lh - mh * q, 2) *
                (1 + lh - mh * q) *
                (deltap + curvature_aver * R + rdeltapp + deltap * shat) -
            2 * epsilon_h_t * (eta - eta_k) * lh * (-1 + lh - mh * q) *
                (lh - mh * q) * (1 + lh - mh * q) * shat *
                std::cos(eta * lh - alpha_0 * mh - eta * mh * q) -
            2 * std::pow(lh, 2) * std::sin(eta) +
            2 * std::pow(lh, 4) * std::sin(eta) +
            4 * lh * mh * q * std::sin(eta) -
            8 * std::pow(lh, 3) * mh * q * std::sin(eta) -
            2 * std::pow(mh, 2) * std::pow(q, 2) * std::sin(eta) +
            12 * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 2) *
                std::sin(eta) -
            8 * lh * std::pow(mh, 3) * std::pow(q, 3) * std::sin(eta) +
            2 * std::pow(mh, 4) * std::pow(q, 4) * std::sin(eta) -
            2 * std::pow(lh, 2) * shat * std::sin(eta) +
            2 * std::pow(lh, 4) * shat * std::sin(eta) +
            4 * lh * mh * q * shat * std::sin(eta) -
            8 * std::pow(lh, 3) * mh * q * shat * std::sin(eta) -
            2 * std::pow(mh, 2) * std::pow(q, 2) * shat * std::sin(eta) +
            12 * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 2) * shat *
                std::sin(eta) -
            8 * lh * std::pow(mh, 3) * std::pow(q, 3) * shat * std::sin(eta) +
            2 * std::pow(mh, 4) * std::pow(q, 4) * shat * std::sin(eta) +
            std::cos(eta) *
                (-2 * (eta - eta_k) * (-1 + lh - mh * q) *
                     std::pow(lh - mh * q, 2) * (1 + lh - mh * q) * shat -
                 (-std::pow(lh, 2) + std::pow(lh, 4) -
                  std::pow(mh, 2) * std::pow(q, 2) +
                  std::pow(mh, 4) * std::pow(q, 4)) *
                     (deltap + rdeltapp + deltap * shat) * std::sin(eta)) -
            deltap * lh * mh * q * std::sin(2 * eta) +
            2 * deltap * std::pow(lh, 3) * mh * q * std::sin(2 * eta) -
            3 * deltap * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 2) *
                std::sin(2 * eta) +
            2 * deltap * lh * std::pow(mh, 3) * std::pow(q, 3) *
                std::sin(2 * eta) -
            lh * mh * q * rdeltapp * std::sin(2 * eta) +
            2 * std::pow(lh, 3) * mh * q * rdeltapp * std::sin(2 * eta) -
            3 * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 2) * rdeltapp *
                std::sin(2 * eta) +
            2 * lh * std::pow(mh, 3) * std::pow(q, 3) * rdeltapp *
                std::sin(2 * eta) -
            deltap * lh * mh * q * shat * std::sin(2 * eta) +
            2 * deltap * std::pow(lh, 3) * mh * q * shat * std::sin(2 * eta) -
            3 * deltap * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 2) *
                shat * std::sin(2 * eta) +
            2 * deltap * lh * std::pow(mh, 3) * std::pow(q, 3) * shat *
                std::sin(2 * eta) +
            deltap * epsilon_h_t * std::pow(lh, 3) *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) -
            deltap * epsilon_h_t * std::pow(lh, 4) *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) -
            2 * deltap * epsilon_h_t * std::pow(lh, 2) * mh * q *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            3 * deltap * epsilon_h_t * std::pow(lh, 3) * mh * q *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            deltap * epsilon_h_t * lh * std::pow(mh, 2) * std::pow(q, 2) *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) -
            3 * deltap * epsilon_h_t * std::pow(lh, 2) * std::pow(mh, 2) *
                std::pow(q, 2) *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            deltap * epsilon_h_t * lh * std::pow(mh, 3) * std::pow(q, 3) *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            epsilon_h_t * std::pow(lh, 3) * rdeltapp *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) -
            epsilon_h_t * std::pow(lh, 4) * rdeltapp *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) -
            2 * epsilon_h_t * std::pow(lh, 2) * mh * q * rdeltapp *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            3 * epsilon_h_t * std::pow(lh, 3) * mh * q * rdeltapp *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            epsilon_h_t * lh * std::pow(mh, 2) * std::pow(q, 2) * rdeltapp *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) -
            3 * epsilon_h_t * std::pow(lh, 2) * std::pow(mh, 2) *
                std::pow(q, 2) * rdeltapp *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            epsilon_h_t * lh * std::pow(mh, 3) * std::pow(q, 3) * rdeltapp *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            deltap * epsilon_h_t * std::pow(lh, 3) * shat *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) -
            deltap * epsilon_h_t * std::pow(lh, 4) * shat *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) -
            2 * deltap * epsilon_h_t * std::pow(lh, 2) * mh * q * shat *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            3 * deltap * epsilon_h_t * std::pow(lh, 3) * mh * q * shat *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            deltap * epsilon_h_t * lh * std::pow(mh, 2) * std::pow(q, 2) *
                shat * std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) -
            3 * deltap * epsilon_h_t * std::pow(lh, 2) * std::pow(mh, 2) *
                std::pow(q, 2) * shat *
                std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) +
            deltap * epsilon_h_t * lh * std::pow(mh, 3) * std::pow(q, 3) *
                shat * std::sin(eta + eta * lh - alpha_0 * mh - eta * mh * q) -
            deltap * epsilon_h_t * std::pow(lh, 3) *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            deltap * epsilon_h_t * std::pow(lh, 4) *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) +
            2 * deltap * epsilon_h_t * std::pow(lh, 2) * mh * q *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) +
            3 * deltap * epsilon_h_t * std::pow(lh, 3) * mh * q *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            deltap * epsilon_h_t * lh * std::pow(mh, 2) * std::pow(q, 2) *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            3 * deltap * epsilon_h_t * std::pow(lh, 2) * std::pow(mh, 2) *
                std::pow(q, 2) *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) +
            deltap * epsilon_h_t * lh * std::pow(mh, 3) * std::pow(q, 3) *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            epsilon_h_t * std::pow(lh, 3) * rdeltapp *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            epsilon_h_t * std::pow(lh, 4) * rdeltapp *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) +
            2 * epsilon_h_t * std::pow(lh, 2) * mh * q * rdeltapp *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) +
            3 * epsilon_h_t * std::pow(lh, 3) * mh * q * rdeltapp *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            epsilon_h_t * lh * std::pow(mh, 2) * std::pow(q, 2) * rdeltapp *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            3 * epsilon_h_t * std::pow(lh, 2) * std::pow(mh, 2) *
                std::pow(q, 2) * rdeltapp *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) +
            epsilon_h_t * lh * std::pow(mh, 3) * std::pow(q, 3) * rdeltapp *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            deltap * epsilon_h_t * std::pow(lh, 3) * shat *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            deltap * epsilon_h_t * std::pow(lh, 4) * shat *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) +
            2 * deltap * epsilon_h_t * std::pow(lh, 2) * mh * q * shat *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) +
            3 * deltap * epsilon_h_t * std::pow(lh, 3) * mh * q * shat *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            deltap * epsilon_h_t * lh * std::pow(mh, 2) * std::pow(q, 2) *
                shat * std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            3 * deltap * epsilon_h_t * std::pow(lh, 2) * std::pow(mh, 2) *
                std::pow(q, 2) * shat *
                std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) +
            deltap * epsilon_h_t * lh * std::pow(mh, 3) * std::pow(q, 3) *
                shat * std::sin(eta - eta * lh + alpha_0 * mh + eta * mh * q) -
            2 * epsilon_h_t * lh * (-1 + lh - mh * q) * (1 + lh - mh * q) *
                (lh - mh * q + shat) *
                std::sin(alpha_0 * mh - eta * (lh - mh * q))) /
           (2. * (-1 + lh - mh * q) * std::pow(lh - mh * q, 2) *
            (1 + lh - mh * q));
}

double Cylinder::g_integration_f(double eta) const {
    return eta;
}
