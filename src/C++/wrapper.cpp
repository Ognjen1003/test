#include "SolveResult.h"
#include "ComponentData.h"
#include "PengRobinson.h"
#include "RachfordRice.h"
#include <omp.h>

std::vector<SolveResult> mainFromPython(const std::vector<Component>& components,
    const std::vector<double>& temperatures,
    const std::vector<double>& pressures) {

    std::vector<SolveResult> results;
    int numTasks = temperatures.size() * pressures.size();
    results.resize(numTasks);

#pragma omp parallel for 
    for (int i = 0; i < temperatures.size(); ++i) {
        for (int j = 0; j < pressures.size(); ++j) {
            double T = temperatures[i];
            double P = pressures[j];
            PengRobinsonEOS eos(components, T, P);
            auto result = RachfordRiceNS::Solver::solve(eos, components, T, P);
            double V = std::get<0>(result);
            int iterations = std::get<3>(result);

            int index = i * pressures.size() + j;
            results[index] = { T, P, V, iterations };
        }
    }

    return results;
}