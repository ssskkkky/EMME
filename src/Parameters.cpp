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
    return -((alpha * eta) / 2.0) - shat * eta * std::cos(eta) + std::sin(eta) +
           shat * std::sin(eta) + 0.25 * alpha * std::sin(2.0 * eta);
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

    return 1.0 / lambda_f_tau_term * bessel_0; // exp_term is throwed in the main part to handle with float overflow
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


    // exp_term is throwed in the main part to handle with float overflow
    std::complex<double> term1 = -bessel_1 *
                                 std::sqrt(bi(eta) * bi(eta_p)) /
                                 pow(lambda_f_tau_term, 3);
    std::complex<double> term2 = bessel_0 * (bi(eta) + bi(eta_p)) /
                                 (2.0 * pow(lambda_f_tau_term, 3));
    std::complex<double> term3 =
        -bessel_0 / pow(lambda_f_tau_term, 2);
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

	auto safe_exp = [](auto var){
	  if (std::real(var)<-40){
	    return std::complex<double> {0.};
	  }
	  else
	  {
	    return std::exp(var);
          }
	};

        return std::pow((q * R) / (taut * vt) * (eta - eta_p), m) *
               safe_exp(
                   -0.5 * std::pow((q * R * (eta - eta_p)) / (vt * taut), 2) +
                   1.i * (beta_1(eta, eta_p) / 2.0 *
                          (((-q * R * (eta - eta_p)) / (vt * taut)))) +
                   std::complex<double>(0.0, 1.0) * taut * omega +
                   -(bi(eta) + bi(eta_p)) /
                       (2.0 * (1.0 + 0.5 * std::complex<double>(0.0, 1.0) *
                                         (taut * vt) / (q * R * (eta - eta_p)) *
                                         beta_1(eta, eta_p)))) /
               taut * (term1 + term2) * jacob;
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
