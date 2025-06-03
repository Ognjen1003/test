import numpy as np
from scipy.optimize import minimize


def gibbs_energy(z, components, eos_model, phase='vapor'):
    z = np.array(z)
    if np.any(z <= 0) or np.abs(np.sum(z) - 1) > 1e-6:
        return 1e10  #problem, vidjet jel moze kako drukcije

    phi_vals = eos_model.fugacity_coeff(z, phase)

    if np.any(np.array(phi_vals) <= 0):
        return 1e10  #ista stvar, fugacitet 0 ili -

    g_star = np.sum(z * np.log(phi_vals))
    return g_star

# Tangent-Plane
def tangent_plane_test(z_trial, components, eos_model):
    return gibbs_energy(z_trial, components, eos_model)

# Hessova matrica - kriticna tocka ovako???
def hessian(z, components, eos_model, delta=1e-5):
    n = len(z)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            z_forward_i = np.array(z)
            z_backward_i = np.array(z)
            z_forward_j = np.array(z)
            z_backward_j = np.array(z)

            z_forward_i[i] += delta
            z_backward_i[i] -= delta
            z_forward_j[j] += delta
            z_backward_j[j] -= delta

            if i == j:
                H[i, j] = (gibbs_energy(z_forward_i, components, eos_model) - 2 * gibbs_energy(z, components, eos_model) + gibbs_energy(z_backward_i, components, eos_model)) / (delta ** 2)
            else:
                H[i, j] = (gibbs_energy(z_forward_i + delta * (np.eye(n)[j]), components, eos_model)
                           - gibbs_energy(z_forward_i - delta * (np.eye(n)[j]), components, eos_model)
                           - gibbs_energy(z_backward_i + delta * (np.eye(n)[j]), components, eos_model)
                           + gibbs_energy(z_backward_i - delta * (np.eye(n)[j]), components, eos_model)) / (4 * delta ** 2)
    return H


def critical_point_function(z_trial, components, eos_model):
    H = hessian(z_trial, components, eos_model)
    eigenvalues = np.linalg.eigvals(H)
    return np.max(np.abs(eigenvalues))


def find_critical_point(components, eos_model, n_components):
    initial_guess = np.array([1.0 / n_components] * n_components)

    bounds = [(1e-6, 1.0 - 1e-6)] * n_components
    cons = ({'type': 'eq', 'fun': lambda z: np.sum(z) - 1})

    result = minimize(tangent_plane_test, initial_guess, args=(components, eos_model), bounds=bounds, constraints=cons)

    if result.success:
        z_eq = result.x
        print("Nađena ravnotežna točka:")
        for i, comp in enumerate(components):
            print(f"{comp['name']}: {z_eq[i]:.4f}")
        print(f"Gibbsova energija u toj točki: g* = {gibbs_energy(z_eq, components, eos_model):.4f}")
    else:
        print("Optimizacija nije uspjela.")

    critical_result = minimize(critical_point_function, initial_guess, args=(components, eos_model), bounds=bounds, constraints=cons)

    if critical_result.success:
        z_crit = critical_result.x
        print("Kritična točka smjese:")
        for i, comp in enumerate(components):
            print(f"{comp['name']}: {z_crit[i]:.4f}")
    else:
        print("Kritična točka nije pronađena.")


# import EOS
# eos_model = EOS(components, T, P)
# find_critical_point(data_oxyfuel_comp3["components"], eos_model, len(data_oxyfuel_comp3["components"]))
