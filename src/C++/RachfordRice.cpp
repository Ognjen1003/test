#include "RachfordRice.h"
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <omp.h>
#include "calculations.h"
#include "EOSBase.h"

namespace RachfordRiceNS {

	double rachfordRiceFunction(double V, const std::vector<double>& z, const std::vector<double>& K) {
		const double eps = 1e-12;  // tolerancija da izbjegnemo dijeljenje s nulom
		double sum = 0.0;
		size_t n = z.size();

		for (size_t i = 0; i < n; ++i) {
			double denom = 1.0 + V * (K[i] - 1.0);
			if (std::abs(denom) > eps) {
				sum += z[i] * (K[i] - 1.0) / denom;
			}
		}

		return sum;
	}


	 std::tuple<double, std::vector<double>, std::vector<double>, int>
		Solver::solve(EOSBase& eos, const std::vector<Component>& components, double T, double P,
			int max_iter, double tol) {

		std::vector<double> z, K;
		for (const auto& comp : components) {
			z.push_back(comp.fraction);
			K.push_back(Calculations::wilson_K(comp, T, P));
		}

		std::tuple<double, std::vector<double>, std::vector<double>, int> last_valid_solution;

		for (int iteration = 0; iteration < max_iter; ++iteration) {
			double K_min = *std::min_element(K.begin(), K.end());
			double K_max = *std::max_element(K.begin(), K.end());

			double V_min = std::max(0.0, 1.0 / (1.0 - K_max));
			double V_max = std::min(1.0, 1.0 / (1.0 - K_min));

			double f_min = rachfordRiceFunction(V_min, z, K);
			double f_max = rachfordRiceFunction(V_max, z, K);

			if (V_min > V_max || f_min * f_max > 0) {
				std::vector<double> y(z.size());
				#pragma omp parallel for
				for (int i = 0; i < z.size(); ++i) y[i] = K[i] * z[i];
				return std::make_tuple(-1.0, z, y, iteration + 1);
			}

			// Fsolve dio - numerička iteracija (Newton-Secant metoda)
			double V = 0.5;
			for (int i = 0; i < max_iter; ++i) {
				double fV = rachfordRiceFunction(V, z, K);
				double dV = 1e-5;
				double fV_dV = rachfordRiceFunction(V + dV, z, K);
				double df = (fV_dV - fV) / dV;
				if (std::abs(df) < 1e-10) break;
				double V_new = V - fV / df;
				if (std::abs(V_new - V) < tol) {
					V = V_new;
					break;
				}
				V = V_new;
			}

			if (!(0.0 <= V && V <= 1.0)) {
				std::vector<double> y(z.size());
				#pragma omp parallel for
				for (int i = 0; i < z.size(); ++i) y[i] = K[i] * z[i];
				return std::make_tuple(-1.0, z, y, iteration + 1);
			}

			std::vector<double> x(z.size()), y(z.size());
			#pragma omp parallel for
			for (int i = 0; i < z.size(); ++i) {
				x[i] = z[i] / (1 + V * (K[i] - 1));
				y[i] = K[i] * x[i];
			}

			//eos.get_Z_factors(x, y);
			std::vector<double> phi_v = eos.fugacity_coeff(y, "vapor");
			std::vector<double> phi_l = eos.fugacity_coeff(x, "liquid");

			std::vector<double> new_K(z.size());
			#pragma omp parallel for
			for (int i = 0; i < z.size(); ++i) {
				new_K[i] = phi_l[i] / phi_v[i];
			}

			last_valid_solution = std::make_tuple(V, x, y, iteration + 1);

			bool K_converged = true;
			#pragma omp parallel for shared(K_converged)
			for (int i = 0; i < K.size(); ++i) {
				if (std::abs((new_K[i] - K[i]) / K[i]) >= tol) {
					#pragma omp critical
					K_converged = false;
					break;
				}
			}

			if (K_converged) {
				return std::make_tuple(V, x, y, iteration + 1);
			}

			K = new_K;
		}

		return last_valid_solution;
	}


} // namespace RachfordRice
