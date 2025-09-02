from src.Endpoints.EOSModul import perform_eos_calculation
from Models.Component import Component
from matplotlib.colors import ListedColormap, BoundaryNorm
from data.testData import ComponentData
from src.EnumsClasses import SolveMethod, EOSType
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
import sys
import os

bin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
sys.path.insert(0, bin_path)  # stavim ga kao prvi prioritet

print("Current dir:", os.getcwd())
print("sys.path:", sys.path)

import eos_cpp 

#u bin imas pyd file, to je dll sa zvanje iz cpp
cplusplus = True

def check_total_fraction(data, label):
    total = sum(comp["fraction"] for comp in data["components"])
    if abs(total - 1.0) < 1e-6:
        print(f"{label}: OK (sum = {total:.6f})")
    else:
        print(f"{label}: NOT OK (sum = {total:.6f})")

def convert_to_cpp_components(py_components):
    cpp_list = []
    for c in py_components:
        cpp_list.append(
            eos_cpp.Component(
                c.name,
                c.formula,
                c.Mw,
                c.Tc,
                c.Pc,
                c.omega,
                c.fraction,
                c.CpA,
                c.CpB,
                c.CpC,
                c.CpD
            )
        )
    return cpp_list     


#check_total_fraction(ComponentData.data_oxyfuel_comp1, "Oxyfuel Comp 1")
#check_total_fraction(ComponentData.data_oxyfuel_comp2, "Oxyfuel Comp 2")
check_total_fraction(ComponentData.data_example_flow_loop, "data_example_flow_loop_4step_30_moles_of_CO2")

temperatures = np.arange(170, 260, 1)  
pressures = np.arange(1, 60, 1)      
results = pd.DataFrame(index=pressures, columns=temperatures)
resultsIteration = pd.DataFrame(index=pressures, columns=temperatures)


components = []

for comp in ComponentData.data_Domagoj_3_ujutro["components"]:
    component = Component(
        name=comp["name"],
        formula=comp["formula"],
        Mw=comp["Mw"],
        Tc=comp["Tc"],
        Pc=comp["Pc"],
        omega=comp["omega"],
        fraction=comp["fraction"],
        CpA=comp["CpA"],
        CpB=comp["CpB"],
        CpC=comp["CpC"],
        CpD=comp["CpD"]
    )
    components.append(component)




start_time = time.time() 

if cplusplus:
    temp_list = temperatures.tolist()
    press_list = pressures.tolist()
    cpp_components = convert_to_cpp_components(components)
    cpp_results = eos_cpp.mainFromPython(cpp_components, temp_list, press_list)

    for res in cpp_results:
        T = res.T
        P = res.P
        V = res.V
        iterations = res.iterations

        if P in results.index and T in results.columns:
            if V == -1:
                V = np.nan 
            results.at[P, T] = V
            resultsIteration.at[P, T] = iterations

else:
    for Tt in temperatures:
        for Pp in pressures:
            result = perform_eos_calculation(
                components,
                Tt,        
                Pp,
                EOSType.PR,            
                SolveMethod.FSOLVE,
                False
                ) 
            if result["V"] == -1.0:
                results.at[Pp, Tt] = np.nan 
                resultsIteration.at[Pp, Tt] = result["iteration"]
            else:
                results.at[Pp, Tt] = result["V"]
                resultsIteration.at[Pp, Tt] = result["iteration"]
                #print(f"{Tt} ---- {Pp}")

end_time = time.time()  # Kraj mjerenja
elapsed_time = end_time - start_time

#print(results)
#print(resultsIteration)
#results.to_csv("results3.xlsx")
#resultsIteration.to_csv("resultsIT.xlsx")

results_float = results.astype(float)


custom_cmap = sns.color_palette("viridis", as_cmap=True)
data_masked = results_float.mask((results_float == 0) | (results_float == 1), np.nan)


cmap_combined = ListedColormap(["black", "white"])
bounds = [0, 0.5, 1.5]
norm = BoundaryNorm(bounds, cmap_combined.N)

plt.figure(figsize=(14, 8))

ax = sns.heatmap(
    data_masked,
    cmap=custom_cmap,
    cbar=True,
    cbar_kws={'label': '[V]'},
    xticklabels=10,
    yticklabels=False 
)

# ➕ Postavi y tickove ručno na svaki 10 bar (npr. 30, 40, ..., 90)
pressures_rounded = pressures.round(1)  # ako imaš 30.0, 30.1, ...
tick_values = list(range(int(pressures[0]), int(pressures[-1])+1, 10))  # 30, 40, ...
tick_indices = [np.abs(pressures_rounded - val).argmin() for val in tick_values]  # indeksi najbliži traženim vrijednostima
ax.set_yticks(tick_indices)
ax.set_yticklabels(tick_values)

# Isto možeš napraviti i za x-tickove ako želiš cijele temperature
temperatures_rounded = temperatures.round(1)
x_tick_values = list(range(int(temperatures[0]), int(temperatures[-1])+1, 5))
x_tick_indices = [np.abs(temperatures_rounded - val).argmin() for val in x_tick_values]
ax.set_xticks(x_tick_indices)
ax.set_xticklabels(x_tick_values, rotation=45)

# Ostale opcije za stil
plt.gcf().patch.set_facecolor('#f0f0f0')
plt.gca().set_facecolor('#f0f0f0')
plt.gca().invert_yaxis()

plt.xlabel("Temperatura [K]")
plt.ylabel("Tlak [bar]")
plt.title("Flow-loop 4th step CO2 jos 30 moles na to")

print(f"Vrijeme : {elapsed_time:.5f} sek")

plt.tight_layout()
plt.show()




