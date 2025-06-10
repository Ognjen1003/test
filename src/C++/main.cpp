#include <vector>
#include <tuple>
#include <iostream>
#include "RachfordRice.h"
#include "PengRobinson.h"
#include "ComponentData.h"

using namespace RachfordRiceNS;


struct SolveResult {
    double T;
    double P;
    double V;
    int iterations;
    double x0;  // npr. prva tekuća kompozicija
    double y0;  // npr. prva parna kompozicija
};

int main() {

    std::vector<double> temperatures = { 280.0, 300.0, 320.0 };
    std::vector<double> pressures = { 10.0, 20.0, 30.0 };

    std::vector<Component> components = ComponentData::getComponents();
    std::vector<SolveResult> results;

    for (double T : temperatures) {
        for (double P : pressures) {
            PengRobinsonEOS eos(components, T, P);
            std::tuple<double, std::vector<double>, std::vector<double>, int> result = Solver::solve(eos, components, T, P);

            double V = std::get<0>(result);
            std::vector<double> x = std::get<1>(result);
            std::vector<double> y = std::get<2>(result);
            int iterations = std::get<3>(result);

            SolveResult r;
            r.T = T;
            r.P = P;
            r.V = V;
            r.iterations = iterations;
            r.x0 = x.empty() ? 0.0 : x[0];
            r.y0 = y.empty() ? 0.0 : y[0];
            results.push_back(r);
        }
    }

    // Pretpostavka: izvoz rezultata u datoteku / serializaciju ide ovdje
    for (const auto& res : results) {
        std::cout << res.T << "," << res.P << "," << res.V << "," << res.iterations << "," << res.x0 << "," << res.y0 << "\n";
    }

    return 0;
}
