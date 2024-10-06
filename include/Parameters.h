#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <array>
#include <complex>

#include "JsonParser.h"

// Structure to hold simulation parameters
struct Parameters {
   public:
    static const Parameters& generate(const util::json::Value&);

    // Member variables
    double q;
    double shat;
    double tau;
    double epsilon_n;
    double eta_i;
    double eta_e;
    double b_theta;
    double beta_e;
    double R;
    double vt;
    double length;
    double theta;
    int npoints;
    int iteration_step_limit;
    double integration_precision;
    double integration_accuracy;
    int integration_iteration_limit;
    int integration_start_points;
    double arc_coeff;
    double alpha;
    double water_bag_weight_vpara;
    double water_bag_weight_vperp;

    // Additional member variables (if needed)
    double omega_s_i;    // Calculated in constructor
    double omega_s_e;    // Calculated in constructor
    double omega_d_bar;  // Calculated in constructor
    bool drift_center_transformation_switch;
    virtual void parameterInit();
    virtual double g_integration_f(double eta) const;
    virtual double beta_1(double eta, double eta_p) const;
    double beta_1_e(double eta, double eta_p) const;
    virtual double bi(double eta) const;

    std::complex<double> lambda_f_tau(double eta,
                                      double eta_p,
                                      std::complex<double> tau) const;
    std::complex<double> integration_lambda_tau(
        double eta,
        double eta_p,
        std::complex<double> tau,
        std::array<std::complex<double>, 5> int_arg) const;
    std::complex<double> integration_lambda_d_tau(
        double eta,
        double eta_p,
        std::complex<double> tau,
        std::array<std::complex<double>, 5> int_arg) const;

    std::complex<double> h_f_tau(std::complex<double> omega,
                                 std::complex<double> tau) const;

    std::complex<double> kappa_f_tau(unsigned int m,
                                     double eta,
                                     double eta_p,
                                     std::complex<double>) const;

    std::complex<double> kappa_f_tau_e(unsigned int m,
                                       double eta,
                                       double eta_p,
                                       std::complex<double>) const;

   protected:
    // Constructor
    Parameters(const util::json::Value&);

    std::array<std::complex<double>, 5> integration_lambda_arg(
        double eta,
        double eta_p,
        std::complex<double> tau) const;
};

struct Stellarator : public Parameters {
    double eta_k;
    int lh;
    int mh;
    double epsilon_h_t;
    double alpha_0;
    double r_over_R;
    double deltap;
    double beta_e_p;
    double rdeltapp;
    double curvature_aver;
    double bi(double eta) const override;
    double g_integration_f(double eta) const override;
    double sigma_f(double eta) const;
    void parameterInit() override;

    friend Parameters;

   protected:
    Stellarator(const util::json::Value&);
};

struct Cylinder : public Parameters {
    Cylinder(const util::json::Value&);

    double shat_coeff;
    double shat_coeff_f(double sv) const;
    double beta_1(double eta, double eta_p) const override;
};

struct CylinderOld : public Parameters {
    CylinderOld(const util::json::Value&);
    double beta_1(double eta, double eta_p) const override;
};

struct TaloyMagneticDrift : public Parameters {
    TaloyMagneticDrift(const util::json::Value&);
    double g_integration_f(double eta) const override;
};

static_assert(std::is_trivially_destructible_v<Parameters>,
              "Paramters should be trivially destructible.");
static_assert(std::is_trivially_destructible_v<Stellarator>,
              "Stellarator should be trivially destructible.");
static_assert(std::is_trivially_destructible_v<Cylinder>,
              "Cylinder should be trivially destructible.");
static_assert(std::is_trivially_destructible_v<CylinderOld>,
              "CylinderOld should be trivially destructible.");
static_assert(std::is_trivially_destructible_v<Cylinder>,
              "TaloyMD should be trivially destructible.");

#endif
