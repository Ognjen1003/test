#ifndef COMPONENT_H
#define COMPONENT_H

#include <string>

class Component {
public:
    std::string name;
    std::string formula;
    double Mw;
    double Tc;
    double Pc;
    double omega;
    double fraction;
    double CpA;
    double CpB;
    double CpC;
    double CpD;

    Component(
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
        double CpD_);
};

#endif // COMPONENT_H
