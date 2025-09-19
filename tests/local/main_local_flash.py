from src.Endpoints.EOSModul import perform_eos_calculation
from Models.Component import Component
from src.Classes import UtilClass
from data.testData import ComponentData
from src.EnumsClasses import SolveMethod, EOSType
import pandas as pd
import numpy as np
import time
import sys
import os

bin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
sys.path.insert(0, bin_path)  # stavim ga kao prvi prioritet

print("Current dir:", os.getcwd())
print("sys.path:", sys.path)

import eos_cpp 

#u bin imas pyd file, to je dll sa zvanje iz cpp
cplusplus = False
display_iterations = False
adjust_display_iterations = False

# funkcije da izgleda urednije
title_primer = "data_oxyfuel_comp1"  # za prikaz vise, nije elementarno
components = ComponentData.oxyfuel_comp1 # podaci koji se actually prikazuju  

UtilClass.check_total_fraction(components, title_primer)

temperatures = np.arange(150, 400, 1)  
pressures = np.arange(1, 105, 1)      
results = pd.DataFrame(index=pressures, columns=temperatures)
resultsIteration = pd.DataFrame(index=pressures, columns=temperatures)

start_time = time.time() 

if cplusplus:
    temp_list = temperatures.tolist()
    press_list = pressures.tolist()
    cpp_components = UtilClass.convert_to_cpp_components(components)
    cpp_results = eos_cpp.mainFromPython(cpp_components, temp_list, press_list)

    for res in cpp_results:
        if res.P in results.index and res.T in results.columns:
            results.at[res.P, res.T] = res.V
            resultsIteration.at[res.P, res.T] = res.iterations

else:
    for Tt in temperatures:
        for Pp in pressures:
            result = perform_eos_calculation(
                components,
                Tt,        
                Pp,
                EOSType.PR,            
                SolveMethod.FSOLVE,
                True
                ) 
            results.at[Pp, Tt] = result["V"]
            resultsIteration.at[Pp, Tt] = result["iteration"]
            #print(f"{Tt} ---- {Pp}")

end_time = time.time()  # Kraj mjerenja
elapsed_time = end_time - start_time
print(f"Vrijeme : {elapsed_time:.5f} sek")

#print(results)
#print(resultsIteration)
#results.to_csv("results3.xlsx")
#resultsIteration.to_csv("resultsIT.xlsx")

if not display_iterations:
    results_display = results.astype(float)
else:
    results_display = resultsIteration.astype(int)

UtilClass.display(temperatures, pressures, results_display, title_primer, adjust_display_iterations)




