import sys
import os
#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
#sys.path.append(os.path.dirname(__file__))
bin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
sys.path.insert(0, bin_path)  # stavimo ga kao prvi prioritet

print("Current dir:", os.getcwd())
print("sys.path:", sys.path)

import eos_cpp 

print("Module loaded:", eos_cpp)

from src.Endpoints.EOSModul import perform_eos_calculation
from src.Classes.Component import Component
from matplotlib.colors import ListedColormap, BoundaryNorm
from tests.testData import ComponentData
from src.EnumsClasses import SolveMethod, EOSType
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import time






print("====================================START C!======================================")



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


check_total_fraction(ComponentData.data_oxyfuel_comp1, "Oxyfuel Comp 1")
check_total_fraction(ComponentData.data_oxyfuel_comp2, "Oxyfuel Comp 2")
check_total_fraction(ComponentData.data_oxyfuel_comp3, "Oxyfuel Comp 3")

temperatures = np.arange(300, 540, 1)  
pressures = np.arange(1, 181, 1)      
results = pd.DataFrame(index=pressures, columns=temperatures)
resultsIteration = pd.DataFrame(index=pressures, columns=temperatures)


components = []

for comp in ComponentData.data["components"]:
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



# Kemijanje, dal ce C++ popusit ovo?
temp_list = temperatures.tolist()
press_list = pressures.tolist()
cpp_components = convert_to_cpp_components(components)


start_time = time.time() 

cpp_results = eos_cpp.mainFromPython(cpp_components, temp_list, press_list)

end_time = time.time()  # Kraj mjerenja
elapsed_time = end_time - start_time
print(f"Vrijeme : {elapsed_time:.5f} sek")

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



print(results)
print(resultsIteration)
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
    yticklabels=10
)


sns.heatmap(
    results_float.where((results_float == 0) | (results_float == 1)),
    cmap=cmap_combined,
    norm=norm,
    cbar=False,
    xticklabels=10,
    yticklabels=10
)


plt.gcf().patch.set_facecolor('#f0f0f0')
plt.gca().set_facecolor('#f0f0f0')
plt.gca().invert_yaxis()


plt.xlabel("Temperatura [K]")
plt.ylabel("Tlak [bar]")
plt.title("C generated")


#black_patch = mpatches.Patch(color='black', label='1 faza')
#white_patch = mpatches.Patch(color='white', label='1 faza')
#plt.legend(handles=[black_patch, white_patch], loc='upper left', title='Jednofazna podrucja')

plt.tight_layout()
plt.show()


