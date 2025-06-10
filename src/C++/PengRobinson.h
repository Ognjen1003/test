#pragma once
#ifndef PENG_ROBINSON_EOS_H
#define PENG_ROBINSON_EOS_H

#include "EOSBase.h"

class PengRobinsonEOS : public EOSBase {
public:
    PengRobinsonEOS(const std::vector<Component>& components, double T, double P);

    double alpha(const Component& comp) const override;
    double a_formula(const Component& comp) const override;
    double b_formula(const Component& comp) const override;
    std::vector<double> get_coefficients(double A, double B) const override;
};

#endif // PENG_ROBINSON_EOS_H


