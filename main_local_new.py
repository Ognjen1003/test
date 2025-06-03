from Endpoints.EOSModul import perform_eos_calculation
from Classes.Component import ComponentNew
from Classes.EOS.EOSUtil import Calculations
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.patches as mpatches
import seaborn as sns
import sys
from CoolProp.CoolProp import PropsSI
import time

start_time = time.time() 

def check_total_fraction(data, label):
    total = sum(comp["fraction"] for comp in data["components"])
    if abs(total - 1.0) < 1e-6:
        print(f"{label}: OK (sum = {total:.6f})")
    else:
        print(f"{label}: NOT OK (sum = {total:.6f})")


data = {
    "components": [
        {
            "name": "Nitrogen", "formula": "N2", "Mw": 28.013, "Tc": 126.2, "Pc": 33.944, "omega": 0.04, "fraction": 0.0,
            "CpA": 28.90, "CpB": -0.0001571, "CpC": 8.081e-07, "CpD": 0.0
        },
        {
            "name": "Carbon Dioxide", "formula": "CO2", "Mw": 44.010, "Tc": 304.7, "Pc": 73.866, "omega": 0.225, "fraction": 0.02,
            "CpA": 22.26, "CpB": 5.981e-3, "CpC": -3.501e-6, "CpD": 7469.0
        },
        {
            "name": "Methane", "formula": "C1", "Mw": 16.043, "Tc": 190.6, "Pc": 46.042, "omega": 0.013, "fraction": 0.443,
            "CpA": 19.89, "CpB": 5.024e-2, "CpC": 1.269e-5, "CpD": -11.01e-9
        },
        {
            "name": "Ethane", "formula": "C2", "Mw": 30.070, "Tc": 305.43, "Pc": 48.839, "omega": 0.0986, "fraction": 0.1407,
            "CpA": 21.13, "CpB": 78.60e-3, "CpC": -11.85e-6, "CpD": 18.99e-9
        },
        {
            "name": "Propane", "formula": "C3", "Mw": 44.097, "Tc": 369.8, "Pc": 42.455, "omega": 0.1524, "fraction": 0.072,
            "CpA": 14.67, "CpB": 75.83e-3, "CpC": -25.60e-6, "CpD": 29890.0
        },
        {
            "name": "i-Butane", "formula": "IC4", "Mw": 58.124, "Tc": 408.1, "Pc": 36.477, "omega": 0.1848, "fraction": 0.038,
            "CpA": 19.26, "CpB": 98.72e-3, "CpC": -37.20e-6, "CpD": 37400.0
        },
        {
            "name": "n-Butane", "formula": "NC4", "Mw": 58.124, "Tc": 425.2, "Pc": 37.966, "omega": 0.201, "fraction": 0.026,
            "CpA": 24.86, "CpB": 133.4e-3, "CpC": -30.1e-6, "CpD": 34.1e-9
        },
        {
            "name": "i-Pentane", "formula": "IC5", "Mw": 72.151, "Tc": 460.4, "Pc": 33.893, "omega": 0.227, "fraction": 0.041,
            "CpA": 23.54, "CpB": 133.7e-3, "CpC": -45.00e-6, "CpD": 47.90e3
        },
        {
            "name": "n-Pentane", "formula": "NC5", "Mw": 72.151, "Tc": 469.6, "Pc": 33.701, "omega": 0.251, "fraction": 0.057,
            "CpA": 23.71, "CpB": 147.4e-3, "CpC": -50.30e-6, "CpD": 52.80e3
        },
        {
            "name": "Hexane", "formula": "C6", "Mw": 84.000, "Tc": 507.5, "Pc": 30.104, "omega": 0.299, "fraction": 0.032,
            "CpA": 25.73, "CpB": 171.8e-3, "CpC": -60.10e-6, "CpD": 61.70e3
        },
        {
            "name": "Heptanes+", "formula": "C7+", "Mw": 156.52, "Tc": 654.31, "Pc": 22.363, "omega": 0.5072, "fraction": 0.1303,
            "CpA": 29.00, "CpB": 220.0e-3, "CpC": -90.00e-6, "CpD": 80.00e3
        }
    ]
}


data_oxyfuel_comp1 = {
    "components": [
        {"name": "CO2", "formula": "CO2", "Mw": 44.01, "Tc": 304.7, "Pc": 73.866, "omega": 0.225, "fraction": 0.85, "CpA": 22.26, "CpB": 5.981e-3, "CpC": -3.501e-6, "CpD": 7469.0},
        {"name": "O2", "formula": "O2", "Mw": 32.00, "Tc": 154.6, "Pc": 50.43, "omega": 0.0222, "fraction": 0.047, "CpA": 29.10, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "N2", "formula": "N2", "Mw": 28.013, "Tc": 126.2, "Pc": 33.944, "omega": 0.04, "fraction": 0.058, "CpA": 28.90, "CpB": -0.0001571, "CpC": 8.081e-7, "CpD": 0.0},
        {"name": "Ar", "formula": "Ar", "Mw": 39.948, "Tc": 150.8, "Pc": 48.63, "omega": -0.004, "fraction": 0.0447, "CpA": 20.85, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "H2O", "formula": "H2O", "Mw": 18.015, "Tc": 647.1, "Pc": 220.64, "omega": 0.344, "fraction": 0.0001, "CpA": 32.24, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "NO", "formula": "NO", "Mw": 30.006, "Tc": 180.0, "Pc": 64.0, "omega": 0.584, "fraction": 0.0001, "CpA": 29.10, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "SO2", "formula": "SO2", "Mw": 64.066, "Tc": 430.8, "Pc": 77.8, "omega": 0.251, "fraction": 0.00005, "CpA": 40.0, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "SO3", "formula": "SO3", "Mw": 80.063, "Tc": 490.0, "Pc": 83.0, "omega": 0.315, "fraction": 0.00005, "CpA": 50.0, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "CO", "formula": "CO", "Mw": 28.010, "Tc": 132.9, "Pc": 34.94, "omega": 0.045, "fraction": 0.00005, "CpA": 29.14, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "H2S", "formula": "H2S+COS", "Mw": 34.08, "Tc": 373.2, "Pc": 88.0, "omega": 0.1, "fraction": 0.00005, "CpA": 34.0, "CpB": 0, "CpC": 0, "CpD": 0}
    ]
}


data_oxyfuel_comp2 = {
    "components": [
        {"name": "CO2", "formula": "CO2", "Mw": 44.01, "Tc": 304.7, "Pc": 73.866, "omega": 0.225, "fraction": 0.98, "CpA": 22.26, "CpB": 5.981e-3, "CpC": -3.501e-6, "CpD": 7469.0},
        {"name": "O2", "formula": "O2", "Mw": 32.00, "Tc": 154.6, "Pc": 50.43, "omega": 0.0222, "fraction": 0.0067, "CpA": 29.10, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "N2", "formula": "N2", "Mw": 28.013, "Tc": 126.2, "Pc": 33.944, "omega": 0.04, "fraction": 0.0071, "CpA": 28.90, "CpB": -0.0001571, "CpC": 8.081e-7, "CpD": 0.0},
        {"name": "Ar", "formula": "Ar", "Mw": 39.948, "Tc": 150.8, "Pc": 48.63, "omega": -0.004, "fraction": 0.0059, "CpA": 20.85, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "H2O", "formula": "H2O", "Mw": 18.015, "Tc": 647.1, "Pc": 220.64, "omega": 0.344, "fraction": 0.0001, "CpA": 32.24, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "NO", "formula": "NO", "Mw": 30.006, "Tc": 180.0, "Pc": 64.0, "omega": 0.584, "fraction": 0.0001, "CpA": 29.10, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "SO2", "formula": "SO2", "Mw": 64.066, "Tc": 430.8, "Pc": 77.8, "omega": 0.251, "fraction": 0.00005, "CpA": 40.0, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "SO3", "formula": "SO3", "Mw": 80.063, "Tc": 490.0, "Pc": 83.0, "omega": 0.315, "fraction": 0.00005, "CpA": 50.0, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "CO", "formula": "CO", "Mw": 28.010, "Tc": 132.9, "Pc": 34.94, "omega": 0.045, "fraction": 0.00005, "CpA": 29.14, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "H2S", "formula": "H2S+COS", "Mw": 34.08, "Tc": 373.2, "Pc": 88.0, "omega": 0.1, "fraction": 0.00005, "CpA": 34.0, "CpB": 0, "CpC": 0, "CpD": 0}
    ]
}


data_oxyfuel_comp3 = {
    "components": [
        {"name": "CO2", "formula": "CO2", "Mw": 44.01, "Tc": 304.7, "Pc": 73.866, "omega": 0.225, "fraction": 0.9994, "CpA": 22.26, "CpB": 5.981e-3, "CpC": -3.501e-6, "CpD": 7469.0},
        {"name": "O2", "formula": "O2", "Mw": 32.00, "Tc": 154.6, "Pc": 50.43, "omega": 0.0222, "fraction": 0.0001, "CpA": 29.10, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "N2", "formula": "N2", "Mw": 28.013, "Tc": 126.2, "Pc": 33.944, "omega": 0.04, "fraction": 0.0001, "CpA": 28.90, "CpB": -0.0001571, "CpC": 8.081e-7, "CpD": 0.0},
        {"name": "Ar", "formula": "Ar", "Mw": 39.948, "Tc": 150.8, "Pc": 48.63, "omega": -0.004, "fraction": 0.0001, "CpA": 20.85, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "H2O", "formula": "H2O", "Mw": 18.015, "Tc": 647.1, "Pc": 220.64, "omega": 0.344, "fraction": 0.0001, "CpA": 32.24, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "NO", "formula": "NO", "Mw": 30.006, "Tc": 180.0, "Pc": 64.0, "omega": 0.584, "fraction": 0.0001, "CpA": 29.10, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "SO2", "formula": "SO2", "Mw": 64.066, "Tc": 430.8, "Pc": 77.8, "omega": 0.251, "fraction": 0.00005, "CpA": 40.0, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "SO3", "formula": "SO3", "Mw": 80.063, "Tc": 490.0, "Pc": 83.0, "omega": 0.315, "fraction": 0.00005, "CpA": 50.0, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "CO", "formula": "CO", "Mw": 28.010, "Tc": 132.9, "Pc": 34.94, "omega": 0.045, "fraction": 0.00005, "CpA": 29.14, "CpB": 0, "CpC": 0, "CpD": 0},
        {"name": "H2S", "formula": "H2S+COS", "Mw": 34.08, "Tc": 373.2, "Pc": 88.0, "omega": 0.1, "fraction": 0.00005, "CpA": 34.0, "CpB": 0, "CpC": 0, "CpD": 0}
    ]
}


check_total_fraction(data_oxyfuel_comp1, "Oxyfuel Comp 1")
check_total_fraction(data_oxyfuel_comp2, "Oxyfuel Comp 2")
check_total_fraction(data_oxyfuel_comp3, "Oxyfuel Comp 3")

temperatures = np.arange(250, 550, 1)  
pressures = np.arange(1, 185, 1)      
results = pd.DataFrame(index=pressures, columns=temperatures)
resultsIteration = pd.DataFrame(index=pressures, columns=temperatures)



components = []

for comp in data["components"]:
    component = ComponentNew(
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


fractions = [comp.fraction for comp in components]


""" result = perform_eos_calculation(
            components,
             288.15,        
             70,            
             fractions,
             "PR"
             ) 
print(f"{result["V"]}")


sys.exit(0) """


for Tt in temperatures:
    for Pp in pressures:
        result = perform_eos_calculation(
            components,
            Tt,        
            Pp,            
            fractions,
            "PR"
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
results.to_csv("results3.xlsx")
resultsIteration.to_csv("resultsIT.xlsx")

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
plt.title("")


#black_patch = mpatches.Patch(color='black', label='1 faza')
#white_patch = mpatches.Patch(color='white', label='1 faza')
#plt.legend(handles=[black_patch, white_patch], loc='upper left', title='Jednofazna podrucja')

plt.tight_layout()
plt.show()


end_time = time.time()  # Kraj mjerenja

elapsed_time = end_time - start_time
print(f"Vrijeme : {elapsed_time:.5f} sek")

#print(f"-------{result["V"]}")



""" components = ["CO2", "O2", "N2", "Ar", "H2O", "CO", "SO2", "NO"]  # CoolProp ne zna za "SO3", "H2S+COS" kao smjesu
for comp in components:
    try:
        Tc = PropsSI("TCRIT", "", 0, "", 0, comp)
        Pc = PropsSI("PCRIT", "", 0, "", 0, comp)
        omega = PropsSI("ACENTRIC", "", 0, "", 0, comp)
        Mw = PropsSI("M", "", 0, "", 0, comp)

        print(f"{comp}: Tc = {Tc:.2f} K, Pc = {Pc/1e5:.2f} bar, omega = {omega:.4f}, M = {Mw:.3f} kg/mol")
    except Exception as e:
        print(f"{comp}: Nije dostupno u CoolProp ({e})")

sys.exit(0)     """   




