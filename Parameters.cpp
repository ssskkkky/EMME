#include "Parameters.h"
#include <cmath>
#include <tuple>
#include "functions.h"

Parameters::Parameters(double q_input,
                       double shat_input,
                       double tau_input,
                       double epsilon_n_input,
                       double eta_i_input,
                       double b_theta_input,
                       double R_input,
                       double vt_input,
                       double length_input,
                       double theta_input,
                       int npoints_input,
                       int iteration_step_limit_input)
    : q(q_input),
      shat(shat_input),
      tau(tau_input),
      epsilon_n(epsilon_n_input),
      eta_i(eta_i_input),
      b_theta(b_theta_input),
      R(R_input),
      vt(vt_input),
      length(length_input),
      theta(theta_input),
      npoints(npoints_input),
      iteration_step_limit(iteration_step_limit_input),
      omega_s_i(-(std::sqrt(b_theta) * vt) / (epsilon_n * R)),
      omega_d_bar(2.0 * epsilon_n * omega_s_i) {}

double Parameters::g_integration_f(double eta) const {
    double alpha = 0.0;
    return -((alpha * eta) / 2.0) - shat * eta * std::cos(eta) + std::sin(eta) +
           shat * std::sin(eta) + 0.25 * alpha * std::sin(2.0 * eta);
}

double Parameters::beta_1(double eta, double eta_p) const {
    return (q * R) / vt * (2.0 * epsilon_n * omega_s_i) *
           (g_integration_f(eta) - g_integration_f(eta_p));
}

double Parameters::bi(double eta) const {
    double alpha = 0.0;
    return b_theta * (1.0 + pow(shat * eta - alpha * std::sin(eta), 2));
}
std::complex<double> Parameters::lambda_f_tau(double eta,
                                              double eta_p,
                                              double tau) const {
    return 1.0 + 0.5 * std::complex<double>(0.0, 1.0) * (tau * vt) /
                     (q * R * (eta - eta_p)) * beta_1(eta, eta_p);
}

// std::array<std::complex<double>, 5>
void Parameters::integration_lambda_arg(double eta, double eta_p, double tau) {
    lambda_f_tau_term = lambda_f_tau(eta, eta_p, tau);
    arg = std::sqrt(bi(eta) * bi(eta_p)) / lambda_f_tau_term;
    // std::complex<double> bessel_0 = util::bessel_ic(0, arg);
    // std::complex<double> bessel_1 = util::bessel_ic(1, arg);
    exp_term = std::exp(-(bi(eta) + bi(eta_p)) / (2.0 * lambda_f_tau_term));
    auto binding_bessel = util::bessel_i_helper(arg);
    bessel_0 = binding_bessel[0];
    bessel_1 = binding_bessel[1];
    return;
}

std::complex<double> Parameters::integration_lambda_tau(double eta,
                                                        double eta_p,
                                                        double tau) const {
    return 1.0 / lambda_f_tau_term * bessel_0 * exp_term;
}

std::complex<double> Parameters::integration_lambda_d_tau(double eta,
                                                          double eta_p,
                                                          double tau) const {
    // auto [arg, bessel_0, bessel_1, exp_term, lambda_f_tau_term] =
    //     integration_lambda_arg(eta, eta_p, tau);
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
                                         double tau) const {
    return exp(std::complex<double>(0.0, 1.0) * tau * omega);
}

std::complex<double> Parameters::kappa_f_tau(double eta,
                                             double eta_p,
                                             std::complex<double> omega)

{
    // Define the integrand function
    auto integrand = [&](double taut) {
        // std::complex<double> arg =
        //     std::sqrt(bi(eta) * bi(eta_p)) / lambda_f_taut(eta, eta_p, taut);

        integration_lambda_arg(eta, eta_p, taut);

        std::complex<double> term1 =
            (omega -
             omega_s_i *
                 (1.0 +
                  eta_i * (0.5 * std::pow((q * R * (eta - eta_p)) / (vt * taut),
                                          2) -
                           1.5))) *
            integration_lambda_tau(eta, eta_p, taut);
        std::complex<double> term2 =
            omega_s_i * eta_i * integration_lambda_d_tau(eta, eta_p, taut);
        return std::exp(-0.5 *
                        std::pow((q * R * (eta - eta_p)) / (vt * taut), 2)) /
               taut * (term1 + term2) *
               std::exp(std::complex<double>(0, 1.0) *
                        (beta_1(eta, eta_p) / 2.0 *
                         (((-q * R * (eta - eta_p)) / (vt * taut))))) *
               h_f_tau(omega, taut);
    };

    // Perform numerical integration
    double lower_bound = 0.0;
    double upper_bound =
        10.0 * M_PI / vt;  // Assuming your integration limits are 0 and 4*pi/vt
    double tol = 1e-5;
    int max_iterations =
        10;  // Adjust as needed for your integration accuracy requirements

    auto result = util::integrate(integrand, lower_bound, upper_bound, tol,
                                  max_iterations);

    // std::complex<double> term1 =
    //     (omega -
    //      omega_s_i *
    //          (1.0 +
    //           eta_i * (0.5 * std::pow((q * R * (eta - eta_p)) / (vt * tau),
    //           2) -
    //                    1.5))) *
    //     integration_lambda_tau(eta, eta_p, tau);

    // auto result = integrand(1.0);

    // Replace "integrate" with the appropriate function call from your chosen
    // numerical integration library

    return -std::complex<double>(0, 1.0) * (q * R) /
           (vt * std::sqrt((2.0 * M_PI))) * result;
}
