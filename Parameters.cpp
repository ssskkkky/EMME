#include "Parameters.h"

#include <cmath>
#include <tuple>

#include "functions.h"

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
                       int iteration_step_limit_input)
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
                                              double tau) const {
    return 1.0 + 0.5 * std::complex<double>(0.0, 1.0) * (tau * vt) /
                     (q * R * (eta - eta_p)) * beta_1(eta, eta_p);
}

std::array<std::complex<double>, 5>
Parameters::integration_lambda_arg(double eta, double eta_p, double tau) const {
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
    double tau,
    std::array<std::complex<double>, 5> int_arg) const {
    auto [lambda_f_tau_term, arg, exp_term, bessel_0, bessel_1] = int_arg;

    return 1.0 / lambda_f_tau_term * bessel_0 * exp_term;
}

std::complex<double> Parameters::integration_lambda_d_tau(
    double eta,
    double eta_p,
    double tau,
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
                                         double tau) const {
    return exp(std::complex<double>(0.0, 1.0) * tau * omega);
}

std::complex<double> Parameters::kappa_f_tau(unsigned int m,
                                             double eta,
                                             double eta_p,
                                             std::complex<double> omega) const

{
    // Define the integrand function
    auto integrand = [&](double taut) {
        // std::complex<double> arg =
        //     std::sqrt(bi(eta) * bi(eta_p)) / lambda_f_taut(eta, eta_p, taut);

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
               std::exp(std::complex<double>(0, 1.0) *
                        (beta_1(eta, eta_p) / 2.0 *
                         (((-q * R * (eta - eta_p)) / (vt * taut))))) *
               h_f_tau(omega, taut);
    };

    // Perform numerical integration
    double lower_bound = 0.0;
    double upper_bound = std::numeric_limits<double>::infinity();
    // 10.0 * M_PI /
    // vt;  // Assuming your integration limits are 0 and 10*pi/vt
    double tol = 1e-5;
    int max_iterations =
        5;  // Adjust as needed for your integration accuracy requirements

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
                 iteration_step_limit_input),
      eta_k(eta_k_input),
      lh(lh_input),
      mh(mh_input),
      epsilon_h_t(epsilon_h_t_input),
      alpha_0(alpha_0_input),
      r_over_R(r_over_R_input),
      deltap(-0.25 * alpha),
      beta_e_p(beta_e * (1.0 + eta_e) / (epsilon_n * R)),
      alpha_p(beta_e_p *
              (q * q / epsilon_n * (1. + eta_e + (1. + eta_i) / tau))),
      deltapp(-0.25 * alpha_p),
      curvature_aver(mh / lh * r_over_R / (q * R) * (4.0 - shat) +
                     (-alpha + 2 * shat * deltap + 0) / R) {}

double Stellarator::sigma_f(double eta) const {
    return shat * (eta - eta_k) +
           (deltap * (1 + shat) + deltapp * epsilon_n * R) * std::sin(eta);
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
    alpha_p =
        beta_e_p * (q * q / epsilon_n * (1. + eta_e + (1. + eta_i) / tau));
    deltapp = -0.25 * alpha_p;
    curvature_aver = mh / lh * r_over_R / (q * R) * (4.0 - shat) +
                     (-alpha + 2 * shat * deltap + 0) / R;
};

double Stellarator::g_integration_f(double eta) const {
    return -0.5 *
           (deltap * std::pow(lh, 3) * q * eta -
            deltap * std::pow(lh, 5) * q * eta -
            2 * deltap * std::pow(lh, 2) * mh * std::pow(q, 2) * eta +
            4 * deltap * std::pow(lh, 4) * mh * std::pow(q, 2) * eta +
            deltap * lh * std::pow(mh, 2) * std::pow(q, 3) * eta -
            6 * deltap * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                eta +
            4 * deltap * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) *
                eta -
            deltap * lh * std::pow(mh, 4) * std::pow(q, 5) * eta +
            4 * std::pow(lh, 2) * mh * r_over_R * eta -
            4 * std::pow(lh, 4) * mh * r_over_R * eta -
            8 * lh * std::pow(mh, 2) * q * r_over_R * eta +
            16 * std::pow(lh, 3) * std::pow(mh, 2) * q * r_over_R * eta +
            4 * std::pow(mh, 3) * std::pow(q, 2) * r_over_R * eta -
            24 * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 2) * r_over_R *
                eta +
            16 * lh * std::pow(mh, 4) * std::pow(q, 3) * r_over_R * eta -
            4 * std::pow(mh, 5) * std::pow(q, 4) * r_over_R * eta +
            3 * deltap * std::pow(lh, 3) * q * shat * eta -
            3 * deltap * std::pow(lh, 5) * q * shat * eta -
            6 * deltap * std::pow(lh, 2) * mh * std::pow(q, 2) * shat * eta +
            12 * deltap * std::pow(lh, 4) * mh * std::pow(q, 2) * shat * eta +
            3 * deltap * lh * std::pow(mh, 2) * std::pow(q, 3) * shat * eta -
            18 * deltap * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                shat * eta +
            12 * deltap * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) *
                shat * eta -
            3 * deltap * lh * std::pow(mh, 4) * std::pow(q, 5) * shat * eta -
            std::pow(lh, 2) * mh * r_over_R * shat * eta +
            std::pow(lh, 4) * mh * r_over_R * shat * eta +
            2 * lh * std::pow(mh, 2) * q * r_over_R * shat * eta -
            4 * std::pow(lh, 3) * std::pow(mh, 2) * q * r_over_R * shat * eta -
            std::pow(mh, 3) * std::pow(q, 2) * r_over_R * shat * eta +
            6 * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 2) * r_over_R *
                shat * eta -
            4 * lh * std::pow(mh, 4) * std::pow(q, 3) * r_over_R * shat * eta +
            std::pow(mh, 5) * std::pow(q, 4) * r_over_R * shat * eta -
            std::pow(lh, 3) * q * alpha * eta +
            std::pow(lh, 5) * q * alpha * eta +
            2 * std::pow(lh, 2) * mh * std::pow(q, 2) * alpha * eta -
            4 * std::pow(lh, 4) * mh * std::pow(q, 2) * alpha * eta -
            lh * std::pow(mh, 2) * std::pow(q, 3) * alpha * eta +
            6 * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) * alpha *
                eta -
            4 * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) * alpha *
                eta +
            lh * std::pow(mh, 4) * std::pow(q, 5) * alpha * eta +
            deltapp * std::pow(lh, 3) * q * R * epsilon_n * eta -
            deltapp * std::pow(lh, 5) * q * R * epsilon_n * eta -
            2 * deltapp * std::pow(lh, 2) * mh * std::pow(q, 2) * R *
                epsilon_n * eta +
            4 * deltapp * std::pow(lh, 4) * mh * std::pow(q, 2) * R *
                epsilon_n * eta +
            deltapp * lh * std::pow(mh, 2) * std::pow(q, 3) * R * epsilon_n *
                eta -
            6 * deltapp * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                R * epsilon_n * eta +
            4 * deltapp * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) *
                R * epsilon_n * eta -
            deltapp * lh * std::pow(mh, 4) * std::pow(q, 5) * R * epsilon_n *
                eta +
            2 * std::pow(lh, 2) * q *
                (-lh + std::pow(lh, 3) + mh * q - 3 * std::pow(lh, 2) * mh * q +
                 3 * lh * std::pow(mh, 2) * std::pow(q, 2) -
                 std::pow(mh, 3) * std::pow(q, 3)) *
                shat * epsilon_h_t * (eta - eta_k) *
                std::cos(lh * eta - mh * (alpha_0 + q * eta)) +
            2 * std::pow(lh, 3) * q * std::sin(eta) -
            2 * std::pow(lh, 5) * q * std::sin(eta) -
            4 * std::pow(lh, 2) * mh * std::pow(q, 2) * std::sin(eta) +
            8 * std::pow(lh, 4) * mh * std::pow(q, 2) * std::sin(eta) +
            2 * lh * std::pow(mh, 2) * std::pow(q, 3) * std::sin(eta) -
            12 * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                std::sin(eta) +
            8 * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) *
                std::sin(eta) -
            2 * lh * std::pow(mh, 4) * std::pow(q, 5) * std::sin(eta) +
            2 * std::pow(lh, 3) * q * shat * std::sin(eta) -
            2 * std::pow(lh, 5) * q * shat * std::sin(eta) -
            4 * std::pow(lh, 2) * mh * std::pow(q, 2) * shat * std::sin(eta) +
            8 * std::pow(lh, 4) * mh * std::pow(q, 2) * shat * std::sin(eta) +
            2 * lh * std::pow(mh, 2) * std::pow(q, 3) * shat * std::sin(eta) -
            12 * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) * shat *
                std::sin(eta) +
            8 * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) * shat *
                std::sin(eta) -
            2 * lh * std::pow(mh, 4) * std::pow(q, 5) * shat * std::sin(eta) +
            lh * q * std::cos(eta) *
                (2 * std::pow(lh - mh * q, 2) *
                     (-1 + std::pow(lh, 2) - 2 * lh * mh * q +
                      std::pow(mh, 2) * std::pow(q, 2)) *
                     shat * (eta - eta_k) +
                 (-std::pow(lh, 2) + std::pow(lh, 4) -
                  std::pow(mh, 2) * std::pow(q, 2) +
                  std::pow(mh, 4) * std::pow(q, 4)) *
                     (deltap + deltap * shat + deltapp * R * epsilon_n) *
                     std::sin(eta)) +
            deltap * std::pow(lh, 2) * mh * std::pow(q, 2) * std::sin(2 * eta) -
            2 * deltap * std::pow(lh, 4) * mh * std::pow(q, 2) *
                std::sin(2 * eta) +
            3 * deltap * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                std::sin(2 * eta) -
            2 * deltap * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) *
                std::sin(2 * eta) +
            deltap * std::pow(lh, 2) * mh * std::pow(q, 2) * shat *
                std::sin(2 * eta) -
            2 * deltap * std::pow(lh, 4) * mh * std::pow(q, 2) * shat *
                std::sin(2 * eta) +
            3 * deltap * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                shat * std::sin(2 * eta) -
            2 * deltap * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) *
                shat * std::sin(2 * eta) +
            deltapp * std::pow(lh, 2) * mh * std::pow(q, 2) * R * epsilon_n *
                std::sin(2 * eta) -
            2 * deltapp * std::pow(lh, 4) * mh * std::pow(q, 2) * R *
                epsilon_n * std::sin(2 * eta) +
            3 * deltapp * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                R * epsilon_n * std::sin(2 * eta) -
            2 * deltapp * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) *
                R * epsilon_n * std::sin(2 * eta) -
            2 * std::pow(lh, 3) * q * epsilon_h_t *
                std::sin(mh * alpha_0 - lh * eta + mh * q * eta) +
            2 * std::pow(lh, 5) * q * epsilon_h_t *
                std::sin(mh * alpha_0 - lh * eta + mh * q * eta) +
            2 * std::pow(lh, 2) * mh * std::pow(q, 2) * epsilon_h_t *
                std::sin(mh * alpha_0 - lh * eta + mh * q * eta) -
            6 * std::pow(lh, 4) * mh * std::pow(q, 2) * epsilon_h_t *
                std::sin(mh * alpha_0 - lh * eta + mh * q * eta) +
            6 * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                epsilon_h_t *
                std::sin(mh * alpha_0 - lh * eta + mh * q * eta) -
            2 * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) *
                epsilon_h_t *
                std::sin(mh * alpha_0 - lh * eta + mh * q * eta) -
            2 * std::pow(lh, 2) * q * shat * epsilon_h_t *
                std::sin(mh * alpha_0 - lh * eta + mh * q * eta) +
            2 * std::pow(lh, 4) * q * shat * epsilon_h_t *
                std::sin(mh * alpha_0 - lh * eta + mh * q * eta) -
            4 * std::pow(lh, 3) * mh * std::pow(q, 2) * shat * epsilon_h_t *
                std::sin(mh * alpha_0 - lh * eta + mh * q * eta) +
            2 * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 3) * shat *
                epsilon_h_t *
                std::sin(mh * alpha_0 - lh * eta + mh * q * eta) +
            deltap * std::pow(lh, 4) * q * epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            deltap * std::pow(lh, 5) * q * epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) -
            2 * deltap * std::pow(lh, 3) * mh * std::pow(q, 2) * epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) -
            3 * deltap * std::pow(lh, 4) * mh * std::pow(q, 2) * epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            deltap * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 3) *
                epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            3 * deltap * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) -
            deltap * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) *
                epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            deltap * std::pow(lh, 4) * q * shat * epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            deltap * std::pow(lh, 5) * q * shat * epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) -
            2 * deltap * std::pow(lh, 3) * mh * std::pow(q, 2) * shat *
                epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) -
            3 * deltap * std::pow(lh, 4) * mh * std::pow(q, 2) * shat *
                epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            deltap * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 3) * shat *
                epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            3 * deltap * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                shat * epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) -
            deltap * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) * shat *
                epsilon_h_t *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            deltapp * std::pow(lh, 4) * q * R * epsilon_h_t * epsilon_n *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            deltapp * std::pow(lh, 5) * q * R * epsilon_h_t * epsilon_n *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) -
            2 * deltapp * std::pow(lh, 3) * mh * std::pow(q, 2) * R *
                epsilon_h_t * epsilon_n *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) -
            3 * deltapp * std::pow(lh, 4) * mh * std::pow(q, 2) * R *
                epsilon_h_t * epsilon_n *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            deltapp * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 3) * R *
                epsilon_h_t * epsilon_n *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) +
            3 * deltapp * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                R * epsilon_h_t * epsilon_n *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) -
            deltapp * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) * R *
                epsilon_h_t * epsilon_n *
                std::sin(mh * alpha_0 + eta - lh * eta + mh * q * eta) -
            deltap * std::pow(lh, 4) * q * epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) +
            deltap * std::pow(lh, 5) * q * epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) +
            2 * deltap * std::pow(lh, 3) * mh * std::pow(q, 2) * epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            3 * deltap * std::pow(lh, 4) * mh * std::pow(q, 2) * epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            deltap * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 3) *
                epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) +
            3 * deltap * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            deltap * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) *
                epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            deltap * std::pow(lh, 4) * q * shat * epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) +
            deltap * std::pow(lh, 5) * q * shat * epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) +
            2 * deltap * std::pow(lh, 3) * mh * std::pow(q, 2) * shat *
                epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            3 * deltap * std::pow(lh, 4) * mh * std::pow(q, 2) * shat *
                epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            deltap * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 3) * shat *
                epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) +
            3 * deltap * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                shat * epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            deltap * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) * shat *
                epsilon_h_t *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            deltapp * std::pow(lh, 4) * q * R * epsilon_h_t * epsilon_n *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) +
            deltapp * std::pow(lh, 5) * q * R * epsilon_h_t * epsilon_n *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) +
            2 * deltapp * std::pow(lh, 3) * mh * std::pow(q, 2) * R *
                epsilon_h_t * epsilon_n *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            3 * deltapp * std::pow(lh, 4) * mh * std::pow(q, 2) * R *
                epsilon_h_t * epsilon_n *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            deltapp * std::pow(lh, 2) * std::pow(mh, 2) * std::pow(q, 3) * R *
                epsilon_h_t * epsilon_n *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) +
            3 * deltapp * std::pow(lh, 3) * std::pow(mh, 2) * std::pow(q, 3) *
                R * epsilon_h_t * epsilon_n *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta)) -
            deltapp * std::pow(lh, 2) * std::pow(mh, 3) * std::pow(q, 4) * R *
                epsilon_h_t * epsilon_n *
                std::sin((1 + lh) * eta - mh * (alpha_0 + q * eta))) /
           (lh * q * (-1 + lh - mh * q) * std::pow(lh - mh * q, 2) *
            (1 + lh - mh * q))

        ;
}
