#include "component.h"

Component::Component(
    const std::string& name_,
    const std::string& formula_,
    double Mw_,
    double Tc_,
    double Pc_,
    double omega_,
    double fraction_,
    double CpA_,
    double CpB_,
    double CpC_,
    double CpD_)
    : name(name_), formula(formula_), Mw(Mw_), Tc(Tc_), Pc(Pc_), omega(omega_),
    fraction(fraction_), CpA(CpA_), CpB(CpB_), CpC(CpC_), CpD(CpD_) {
}


