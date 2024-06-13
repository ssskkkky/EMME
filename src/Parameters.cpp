#include "Parameters.h"

#include <cmath>
#include <iostream>
#include <tuple>

#include "functions.h"
using namespace std::literals;
Parameters::Parameters(double q_input,
                       double shat_input,
                       double tau_input,
                       double epsilon_n_input,
                       double eta_i_input,
                       double eta_e_input,
                       double b_theta_input,
                       double beta_e_input,
                       double R_input,
                       double vt_input,
                       double length_input,
                       double theta_input,
                       int npoints_input,
                       int iteration_step_limit_input,
                       double integration_precision_input,
                       double integration_accuracy_input,
                       int integration_iteration_limit_input,
                       int integration_start_points_input,
                       double arc_coeff_input)
    : q(q_input),
      shat(shat_input),
      tau(tau_input),
      epsilon_n(epsilon_n_input),
      eta_i(eta_i_input),
      eta_e(eta_e_input),
      b_theta(b_theta_input),
      beta_e(beta_e_input),
      R(R_input),
      vt(vt_input),
      alpha(q * q * R * beta_e / (epsilon_n * R) *
            ((1 + eta_e) + 1 / tau * (1 + eta_i))),
      length(length_input),
      theta(theta_input),
      npoints(npoints_input),
      iteration_step_limit(iteration_step_limit_input),
      integration_precision(integration_precision_input),
      integration_accuracy(integration_accuracy_input),
      integration_iteration_limit(integration_iteration_limit_input),
      integration_start_points(integration_start_points_input),
      arc_coeff(arc_coeff_input),
      omega_s_i(-(std::sqrt(b_theta) * vt) / (epsilon_n * R)),
      omega_s_e(-tau * omega_s_i),
      omega_d_bar(2.0 * epsilon_n * omega_s_i) {}

void Parameters::parameterInit() {
    alpha = q * q * R * beta_e / (epsilon_n * R) *
            ((1 + eta_e) + 1 / tau * (1 + eta_i));
    omega_s_i = -(std::sqrt(b_theta) * vt) / (epsilon_n * R);
    omega_s_e = -tau * omega_s_i;
    omega_d_bar = 2.0 * epsilon_n * omega_s_i;
};

double Parameters::g_integration_f(double eta) const {
    return -((alpha * eta) / 2.0) - shat * eta * std::cos(eta) + std::sin(eta) +
           shat * std::sin(eta) + 0.25 * alpha * std::sin(2.0 * eta);
}

double Parameters::beta_1(double eta, double eta_p) const {
    return (q * R) / vt * (2.0 * epsilon_n * omega_s_i) *
           (g_integration_f(eta) - g_integration_f(eta_p));
}

double Parameters::beta_1_e(double eta, double eta_p) const {
    return (q * R) / vt * (2.0 * epsilon_n * omega_s_e) *
           (g_integration_f(eta) - g_integration_f(eta_p));
}

double Parameters::bi(double eta) const {
    return b_theta * (1.0 + pow(shat * eta - alpha * std::sin(eta), 2));
}
std::complex<double> Parameters::lambda_f_tau(double eta,
                                              double eta_p,
                                              std::complex<double> tau) const {
    return 1.0 + 0.5 * std::complex<double>(0.0, 1.0) * (tau * vt) /
                     (q * R * (eta - eta_p)) * beta_1(eta, eta_p);
}

std::array<std::complex<double>, 5> Parameters::integration_lambda_arg(
    double eta,
    double eta_p,
    std::complex<double> tau) const {
    std::complex<double> lambda_f_tau_term = lambda_f_tau(eta, eta_p, tau);
    std::complex<double> arg =
        std::sqrt(bi(eta) * bi(eta_p)) / lambda_f_tau_term;
    // std::complex<double> bessel_0 = util::bessel_ic(0, arg);
    // std::complex<double> bessel_1 = util::bessel_ic(1, arg);
    std::complex<double> exp_term =
        std::exp(-(bi(eta) + bi(eta_p)) / (2.0 * lambda_f_tau_term));
    auto binding_bessel = util::bessel_i_helper(arg);
    std::complex<double> bessel_0 = binding_bessel[0];
    std::complex<double> bessel_1 = binding_bessel[1];
    return {lambda_f_tau_term, arg, exp_term, bessel_0, bessel_1};
}

std::complex<double> Parameters::integration_lambda_tau(
    double eta,
    double eta_p,
    std::complex<double> tau,
    std::array<std::complex<double>, 5> int_arg) const {
    auto [lambda_f_tau_term, arg, exp_term, bessel_0, bessel_1] = int_arg;

    return 1.0 / lambda_f_tau_term * bessel_0 * exp_term;
}

std::complex<double> Parameters::integration_lambda_d_tau(
    double eta,
    double eta_p,
    std::complex<double> tau,
    std::array<std::complex<double>, 5> int_arg) const {
    auto [lambda_f_tau_term, arg, exp_term, bessel_0, bessel_1] = int_arg;
    // std::complex<double> arg =
    //     std::sqrt(bi(eta) * bi(eta_p)) / lambda_f_tau(eta, eta_p, tau);
    // std::complex<double> bessel_0 = util::bessel_ic(0, arg);
    // std::complex<double> bessel_1 = util::bessel_ic(1, arg);
    // std::complex<double> exp_term = std::exp(
    //     -(bi(eta) + bi(eta_p)) / (2.0 * lambda_f_tau(eta, eta_p, tau)));

    // std::complex<double> lambda_f_tau_term = lambda_f_tau(eta, eta_p, tau);

    std::complex<double> term1 = -exp_term * bessel_1 *
                                 std::sqrt(bi(eta) * bi(eta_p)) /
                                 pow(lambda_f_tau_term, 3);
    std::complex<double> term2 = exp_term * bessel_0 * (bi(eta) + bi(eta_p)) /
                                 (2.0 * pow(lambda_f_tau_term, 3));
    std::complex<double> term3 =
        -exp_term * bessel_0 / pow(lambda_f_tau_term, 2);
    return term1 + term2 + term3;
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
        // std::complex<double> arg =
        //     std::sqrt(bi(eta) * bi(eta_p)) / lambda_f_taut(eta, eta_p, taut);
        auto omi = -std::copysign(1, omega.real());
        auto taut = arc_coeff * std::atan(taut_transformed) -
                    1.i * omi * taut_transformed;
        auto jacob =
            arc_coeff / (1 + taut_transformed * taut_transformed) - 1.i * omi;

        auto passer = integration_lambda_arg(eta, eta_p, taut);

        std::complex<double> term1 =
            (omega -
             omega_s_i *
                 (1.0 +
                  eta_i * (0.5 * std::pow((q * R * (eta - eta_p)) / (vt * taut),
                                          2) -
                           1.5))) *
            integration_lambda_tau(eta, eta_p, taut, passer);

        std::complex<double> term2 =
            omega_s_i * eta_i *
            integration_lambda_d_tau(eta, eta_p, taut, passer);

        return std::pow((q * R) / (taut * vt) * (eta - eta_p), m) *
               std::exp(-0.5 *
                        std::pow((q * R * (eta - eta_p)) / (vt * taut), 2)) /
               taut * (term1 + term2) *
               std::exp(1.i * (beta_1(eta, eta_p) / 2.0 *
                               (((-q * R * (eta - eta_p)) / (vt * taut))))) *
               h_f_tau(omega, taut) * jacob;
    };

    auto result =
        util::integrate(integrand, integration_precision, integration_accuracy,
                        integration_iteration_limit, integration_start_points);

    // std::complex<double> term1 =
    //     (omega -
    //      omega_s_i *
    //          (1.0 +
    //           eta_i * (0.5 * std::pow((q * R * (eta - eta_p)) / (vt *
    //           tau), 2) -
    //                    1.5))) *
    //     integration_lambda_tau(eta, eta_p, tau);

    // auto result = integrand(1.0);

    // Replace "integrate" with the appropriate function call from your
    // chosen numerical integration library

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

Stellarator::Stellarator(double q_input,
                         double shat_input,
                         double tau_input,
                         double epsilon_n_input,
                         double eta_i_input,
                         double eta_e_input,
                         double b_theta_input,
                         double beta_e_input,
                         double R_input,
                         double vt_input,
                         double length_input,
                         double theta_input,
                         int npoints_input,
                         int iteration_step_limit_input,
                         double integration_precision_input,
                         double integration_accuracy_input,
                         int integration_iteration_limit_input,
                         int integration_start_points_input,
                         double arc_coeff_input,
                         double eta_k_input,
                         int lh_input,
                         int mh_input,
                         double epsilon_h_t_input,
                         double alpha_0_input,
                         double r_over_R_input)
    : Parameters(q_input,
                 shat_input,
                 tau_input,
                 epsilon_n_input,
                 eta_i_input,
                 eta_e_input,
                 b_theta_input,
                 beta_e_input,
                 R_input,
                 vt_input,
                 length_input,
                 theta_input,
                 npoints_input,
                 iteration_step_limit_input,
                 integration_precision_input,
                 integration_accuracy_input,
                 integration_iteration_limit_input,
                 integration_start_points_input,
                 arc_coeff_input),
      eta_k(eta_k_input),
      lh(lh_input),
      mh(mh_input),
      epsilon_h_t(epsilon_h_t_input),
      alpha_0(alpha_0_input),
      r_over_R(r_over_R_input),
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
    omega_d_bar = 2.0 * epsilon_n * omega_s_i;

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

Cylinder::Cylinder(double q_input,
                   double shat_input,
                   double tau_input,
                   double epsilon_n_input,
                   double eta_i_input,
                   double eta_e_input,
                   double b_theta_input,
                   double beta_e_input,
                   double R_input,
                   double vt_input,
                   double length_input,
                   double theta_input,
                   int npoints_input,
                   int iteration_step_limit_input,
                   double integration_precision_input,
                   double integration_accuracy_input,
                   int integration_iteration_limit_input,
                   int integration_start_points_input,
                   double arc_coeff_input)
    : Parameters(q_input,
                 shat_input,
                 tau_input,
                 epsilon_n_input,
                 eta_i_input,
                 eta_e_input,
                 b_theta_input,
                 beta_e_input,
                 R_input,
                 vt_input,
                 length_input,
                 theta_input,
                 npoints_input,
                 iteration_step_limit_input,
                 integration_precision_input,
                 integration_accuracy_input,
                 integration_iteration_limit_input,
                 integration_start_points_input,
                 arc_coeff_input),
      shat_coeff(shat_coeff_f(shat)) {
    std::cout << shat_coeff << std::endl;
}

double Cylinder::shat_coeff_f(double sv) const {
    auto ansx_shat = 3.14973 - 0.374507 / sv + 0.0235059 / (sv * sv);
    return g_integration_f(ansx_shat) / ansx_shat;
}
// double Cylinder::beta_1(double eta) const {
//     return shat * eta;
// }

double Cylinder::beta_1(double eta, double eta_p) const {
    return (q * R) / vt * (2.0 * epsilon_n * omega_s_i) * (eta - eta_p) *
           shat_coeff;
}
