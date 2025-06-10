#include "EOSBase.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

EOSBase::EOSBase(const std::vector<Component>& components_, double T_, double P_)
    : components(components_), T(T_), P(P_) {
}

void EOSBase::calc_parameters() {
    for (const auto& comp : components) {
        double alpha_i = this->alpha(comp);
        a_i.push_back(a_formula(comp) * alpha_i);
        b_i.push_back(b_formula(comp));
    }
}

std::pair<double, double> EOSBase::get_Z_factors(const std::vector<double>& x, const std::vector<double>& y) {
    auto compute_A_B = [&](const std::vector<double>& frac) {
        double a_mix = 0.0;
        double b_mix = 0.0;
        size_t N = frac.size();
        for (size_t i = 0; i < N; ++i) {
            b_mix += frac[i] * b_i[i];
            for (size_t j = 0; j < N; ++j) {
                a_mix += frac[i] * frac[j] * std::sqrt(a_i[i] * a_i[j]);
            }
        }
        double A = a_mix * P / (R * R * T * T);
        double B = b_mix * P / (R * T);
        return std::make_pair(A, B);
        };

    std::pair<double, double> AB_l = compute_A_B(x);
    std::pair<double, double> AB_v = compute_A_B(y);

    Zl = calc_Z_factor(AB_l.first, AB_l.second, "liquid");
    Zv = calc_Z_factor(AB_v.first, AB_v.second, "vapor");

    return { Zl, Zv };
}

double EOSBase::calc_Z_factor(double A, double B, const std::string& phase) {
    std::vector<double> coeffs = get_coefficients(A, B);
    std::vector<double> roots = solve_cubic(coeffs);
    if (roots.empty()) throw std::runtime_error("No real roots found for cubic EOS");
    return (phase == "vapor") ? *std::max_element(roots.begin(), roots.end())
        : *std::min_element(roots.begin(), roots.end());
}

std::vector<double> EOSBase::solve_cubic(const std::vector<double>& coeffs) {
    double a = coeffs[0];
    double b = coeffs[1];
    double c = coeffs[2];
    double d = coeffs[3];

    double f = ((3.0 * c / a) - (b * b / (a * a))) / 3.0;
    double g = ((2.0 * b * b * b / (a * a * a)) - (9.0 * b * c / (a * a)) + (27.0 * d / a)) / 27.0;
    double h = (g * g / 4.0) + (f * f * f / 27.0);

    std::vector<double> roots;
    if (h > 0) {
        double R = -(g / 2.0) + std::sqrt(h);
        double S = std::cbrt(R);
        double T = -(g / 2.0) - std::sqrt(h);
        double U = std::cbrt(T);
        roots.push_back((S + U) - (b / (3.0 * a)));
    }
    else {
        double i = std::sqrt((g * g / 4.0) - h);
        double j = std::cbrt(i);
        double k = std::acos(-(g / (2.0 * i)));
        double M = std::cos(k / 3.0);
        double N = std::sqrt(3.0) * std::sin(k / 3.0);
        double P = -b / (3.0 * a);
        roots.push_back(2.0 * j * M + P);
        roots.push_back(-j * (M + N) + P);
        roots.push_back(-j * (M - N) + P);
    }
    return roots;
}

std::vector<double> EOSBase::fugacity_coeff(const std::vector<double>& x, const std::string& phase) {
    double a_mix = 0.0;
    double b_mix = 0.0;
    size_t N = x.size();

    for (size_t i = 0; i < N; ++i) {
        b_mix += x[i] * b_i[i];
        for (size_t j = 0; j < N; ++j) {
            a_mix += x[i] * x[j] * std::sqrt(a_i[i] * a_i[j]);
        }
    }

    double A = a_mix * P / (R * R * T * T);
    double B = b_mix * P / (R * T);
    double Z = calc_Z_factor(A, B, phase);

    std::vector<double> phi;
    for (size_t i = 0; i < N; ++i) {
        double bi = b_i[i];
        double ai = a_i[i];
        double sum_a = 0.0;
        for (size_t j = 0; j < N; ++j) {
            sum_a += x[j] * std::sqrt(ai * a_i[j]);
        }
        double term1 = bi / b_mix * (Z - 1.0) - std::log(Z - B);
        double term2 = A / (2.0 * std::sqrt(2.0) * B);
        double term3 = (2.0 * sum_a / a_mix) - (bi / b_mix);
        double term4 = std::log((Z + (1.0 + std::sqrt(2.0)) * B) / (Z + (1.0 - std::sqrt(2.0)) * B));
        double ln_phi = term1 - term2 * term3 * term4;
        phi.push_back(std::exp(ln_phi));
    }
    return phi;
}

