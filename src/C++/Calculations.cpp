#include "calculations.h"

double Calculations::wilson_K(const Component& component, double T, double P) {
    double Tr = T / component.Tc;
    return std::exp(std::log(component.Pc / P) + 5.373 * (1 + component.omega) * (1 - 1 / Tr));
}