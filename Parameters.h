#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <complex>
// Structure to hold simulation parameters
struct Parameters {
   public:
    // Constructor
    Parameters(double q_input,
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
               int iteration_step_limit_input);

    // Member variables
    double q;
    double shat;
    double tau;
    double epsilon_n;
    double eta_i;
    double b_theta;
    double R;
    double vt;
    double length;
    double theta;
    int npoints;
    int iteration_step_limit;

    // Additional member variables (if needed)
    double omega_s_i;    // Calculated in constructor
    double omega_d_bar;  // Calculated in constructor
    double g_integration_f(double eta);
    double beta_1(double eta, double eta_p);
    double bi(double eta);

    std::complex<double> lambda_f_tau(double eta, double eta_p, double tau);
    std::complex<double> integration_lambda_tau(double eta,
                                                double eta_p,
                                                double tau);
    std::complex<double> integration_lambda_d_tau(double eta,
                                                  double eta_p,
                                                  double tau);

    std::complex<double> h_f_tau(std::complex<double> omega, double tau);

    std::complex<double> kappa_f_tau(double eta,
                                     double eta_p,
                                     std::complex<double>);

   private:
    // No private member functions needed (assuming this is just a data
    // structure)
};

template <class T>
std::complex<double> integrate(T, double, double, int) {
    return 0;
}  // todo
double bessel_ic(int, std::complex<double>);  // todo
#endif                                        // PARAMETERS_H
