from typing import List, Dict
import src.Classes.EOS.RachfordRice as RR
from src.Classes.EOS import PengRobinsonEOS
from Models.Component import Component
import numpy as np
import matplotlib.pyplot as plt

R = 8.314462618  # J/(mol·K)
V_LOOP = 0.05  # m3

# === DEFINICIJA KOMPONENTI ===
components = [
    Component(name="CarbonDioxide", formula="CO2", Mw=44.01, Tc=304.7, Pc=73.866, omega=0.225, fraction=0.91, CpA=22.26, CpB=5.981e-3, CpC=-3.501e-6, CpD=7469.0),
    Component(name="Nitrogen", formula="N2", Mw=28.013, Tc=126.2, Pc=33.944, omega=0.04, fraction=0.04, CpA=28.90, CpB=-0.0001571, CpC=8.081e-7, CpD=0.0),
    Component(name="Argon", formula="Ar", Mw=39.948, Tc=150.8, Pc=48.63, omega=-0.004, fraction=0.003, CpA=20.85, CpB=0.0, CpC=0.0, CpD=0.0),
    Component(name="Water", formula="H2O", Mw=18.015, Tc=647.1, Pc=220.64, omega=0.344, fraction=0.002, CpA=32.24, CpB=0.0, CpC=0.0, CpD=0.0),
    Component(name="SulfurDioxide", formula="SO2", Mw=64.066, Tc=430.8, Pc=77.8, omega=0.251, fraction=0.045, CpA=40.0, CpB=0.0, CpC=0.0, CpD=0.0)
]

# === CILJNI SASTAV ===
target_composition = {
    "CarbonDioxide": 0.91,
    "Nitrogen": 0.04,
    "SulfurDioxide": 0.045,
    "Argon": 0.003,
    "Water": 0.002
}

def calc_pressure_from_composition(eos_class, components: List, composition: Dict[str, float],
                                    T_K: float, V_m3: float, n_total: float) -> float:
    P_bar_guess = 50.0
    for comp in components:
        comp.fraction = composition.get(comp.name, 0.0)
    result = RR.solve(eos_class, components, T_K, P_bar_guess)
    V_frac, x, y, method, n_iter, Zl, Zv = result
    Z_used = Zv if V_frac > 0.99 else Zl if V_frac < 0.01 else V_frac * Zv + (1 - V_frac) * Zl
    P_Pa = n_total * Z_used * R * T_K / V_m3
    return P_Pa / 1e5  # u bar

# === KORAK 1: POCETNO PUNJENJE CO2 na 1.1 bar i 20°C ===
T0 = 293.15  # 20°C
P0 = 1.1e5   # 1.1 bar
for comp in components:
    comp.fraction = 1.0 if comp.name == "CarbonDioxide" else 0.0
Z0_result = RR.solve(PengRobinsonEOS, components, T0, P0)
_, _, _, _, _, Zl0, Zv0 = Z0_result
Z_CO2 = Zv0
n_CO2_start = P0 * V_LOOP / (Z_CO2 * R * T0)

# === KORAK 2: SEKVENCIJALNO DODAVANJE KOMPONENTI ===
added_components = ["SulfurDioxide", "Argon", "Nitrogen", "Water"]
composition_steps = []
current_mixture = {"CarbonDioxide": n_CO2_start}

for comp_name in added_components:
    target_frac = target_composition[comp_name]
    total_so_far = sum(current_mixture.values())
    n_added = target_frac / (1 - target_composition["CarbonDioxide"]) * total_so_far
    current_mixture[comp_name] = n_added
    composition_steps.append(current_mixture.copy())

# === KORAK 3: NADOPUNA CO2 na takav broj molova da pritisak bude 70 bara pri 60°C ===
T_final = 333.15  # 60°C
P_target_bar = 70.0
from scipy.optimize import fsolve

def co2_to_reach_target_pressure(n_existing: float, mixture: Dict[str, float]) -> float:
    def pressure_residual(n_co2):
        mixture["CarbonDioxide"] = n_co2
        n_total = sum(mixture.values())
        mol_fractions = {k: v / n_total for k, v in mixture.items()}
        p_bar = calc_pressure_from_composition(PengRobinsonEOS, components, mol_fractions, T_final, V_LOOP, n_total)
        return p_bar - P_target_bar

    n_guess = n_existing * 2
    n_solution = fsolve(pressure_residual, n_guess)[0]
    return n_solution

n_CO2_final = co2_to_reach_target_pressure(n_CO2_start, current_mixture)
current_mixture["CarbonDioxide"] = n_CO2_final
composition_steps.append(current_mixture.copy())

# === TEMPERATURE LIST ===
T_list = [T0] * (len(composition_steps) - 1) + [T_final]

# === SIMULACIJA ===
pressures = []
for step, T_K in zip(composition_steps, T_list):
    n_total = sum(step.values())
    mol_fractions = {k: v / n_total for k, v in step.items()}
    p_bar = calc_pressure_from_composition(PengRobinsonEOS, components, mol_fractions, T_K, V_LOOP, n_total)
    pressures.append(p_bar)

# === ISPIS ===
for i, p in enumerate(pressures):
    print(f"p{i} = {p:.2f} bar")

# === GRAF ===
plt.plot([f"p{i}" for i in range(len(pressures))], pressures, marker='o')
plt.ylabel("Tlak [bar]")
plt.title("Sekvenca miješanja: p0 - pN s kompresijom CO2 do 70 bara")
plt.grid(True)
plt.show()
