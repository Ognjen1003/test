#include "PengRobinson.h"
#include <cmath>


PengRobinsonEOS::PengRobinsonEOS(const std::vector<Component>& components, double T, double P) : EOSBase(components, T, P) {
    calc_parameters();  // SIGURNO: alpha() je sada validan virtualni poziv
}

double PengRobinsonEOS::alpha(const Component& comp) const {
    double Tr = T / comp.Tc;
    double omega = comp.omega;
    double m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    return std::pow(1 + m * (1 - std::sqrt(Tr)), 2);
}

double PengRobinsonEOS::a_formula(const Component& comp) const {
    return 0.45724 * std::pow(R, 2) * std::pow(comp.Tc, 2) / comp.Pc;
}

double PengRobinsonEOS::b_formula(const Component& comp) const {
    return 0.07780 * R * comp.Tc / comp.Pc;
}

std::vector<double> PengRobinsonEOS::get_coefficients(double A, double B) const {
    return {
        1.0,
        -(1.0 - B),
        A - 3.0 * B * B - 2.0 * B,
        -(A * B - B * B - B * B * B)
    };
}

