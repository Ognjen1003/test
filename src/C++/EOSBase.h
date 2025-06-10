#pragma once
#ifndef EOS_BASE_H
#define EOS_BASE_H

#include <vector>
#include "component.h"

class EOSBase {
public:
    static constexpr double R = 8.314;

    EOSBase(const std::vector<Component>& components, double T, double P);

    virtual double alpha(const Component& comp) const = 0;
    virtual double a_formula(const Component& comp) const = 0;
    virtual double b_formula(const Component& comp) const = 0;
    virtual std::vector<double> get_coefficients(double A, double B) const = 0;

    double Zl, Zv;

    std::pair<double, double> get_Z_factors(const std::vector<double>& x, const std::vector<double>& y);
    std::vector<double> fugacity_coeff(const std::vector<double>& x, const std::string& phase);

protected:
    std::vector<Component> components;
    std::vector<double> a_i;
    std::vector<double> b_i;
    double T, P;

    void calc_parameters();
    double calc_Z_factor(double A, double B, const std::string& phase);
    std::vector<double> solve_cubic(const std::vector<double>& coeffs);
};

#endif // EOS_BASE_H


