from Endpoints.EOSModul import perform_eos_calculation
from Classes.Component import ComponentNew
from Classes.EOS.EOSUtil import Calculations
import math


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

K = [Calculations.wilson_K(comp, 288, 25) for comp in components]
print(f"[DEBUG] Inicijalni K: {K}")


fractions = [comp.fraction for comp in components]
result = perform_eos_calculation(
    components,
    288,        
    25,            
    fractions,
    "PR"
) 

print(result)




