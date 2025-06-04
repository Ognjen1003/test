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



#import sys
#import os
#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

start_time = time.time() 

def check_total_fraction(data, label):
    total = sum(comp["fraction"] for comp in data["components"])
    if abs(total - 1.0) < 1e-6:
        print(f"{label}: OK (sum = {total:.6f})")
    else:
        print(f"{label}: NOT OK (sum = {total:.6f})")


check_total_fraction(ComponentData.data_oxyfuel_comp1, "Oxyfuel Comp 1")
check_total_fraction(ComponentData.data_oxyfuel_comp2, "Oxyfuel Comp 2")
check_total_fraction(ComponentData.data_oxyfuel_comp3, "Oxyfuel Comp 3")

temperatures = np.arange(300, 540, 2)  
pressures = np.arange(1, 181, 2)      
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
            print(f"{Tt} ---- {Pp}")


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
plt.title("OxyFuel 85% CO2")


#black_patch = mpatches.Patch(color='black', label='1 faza')
#white_patch = mpatches.Patch(color='white', label='1 faza')
#plt.legend(handles=[black_patch, white_patch], loc='upper left', title='Jednofazna podrucja')

plt.tight_layout()
plt.show()


end_time = time.time()  # Kraj mjerenja

elapsed_time = end_time - start_time
print(f"Vrijeme : {elapsed_time:.5f} sek")
