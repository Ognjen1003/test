#ifndef RACHFORD_RICE_H
#define RACHFORD_RICE_H

#include <vector>
#include <tuple>
#include "component.h"
#include "EOSBase.h"

namespace RachfordRiceNS {

    double rachfordRiceFunction(double V, const std::vector<double>& z, const std::vector<double>& K);

    class Solver {
    public:
        static std::tuple<double, std::vector<double>, std::vector<double>, int>
            solve(EOSBase& eos, const std::vector<Component>& components, double T, double P,
                int max_iter = 100, double tol = 1e-6);
    };

} // namespace RachfordRice

#endif // RACHFORD_RICE_H


