from src.Endpoints.EOSModul import perform_eos_calculation
from src.Classes.UtilClass import Util
from data.testData import ComponentData
from data.testDataBIC import ComponentDataBIC
import src.EnumsClasses.MethodsAndTypes as MT
import pandas as pd
import numpy as np
import time, sys, os

bin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
sys.path.insert(0, bin_path)  # stavim ga kao prvi prioritet

print("Current dir:", os.getcwd())
print("sys.path:", sys.path)

import eos_cpp 

#u bin imas pyd file, to je dll sa zvanje iz cpp
cplusplus = False                   # python mora imati istu verziju kao i pyd file (313 npr -> 3.13), pazi da radi sa BIC, BIC True ne radi ovdje
display_iterations = False          # samo iteracije nema veze sa negative flesh nuzno
adjust_display_iterations = False   # prikazuje u razlicitm bojama nizi broj iteracija, bitno i za negativan flash
toggle_phase_detect = False         # umjesto V pare daje sifru izlaza iz Rachford-Rice, sa display_iterations, ostalo zgasi
is_BIC_used = False                  # provjeri koju matricu uopce upotrebljavas


# funkcije da izgleda urednije
title_primer = "oxyfuel_comp1"  # za prikaz vise, nije elementarno
components = ComponentData.data_components # podaci koji se actually prikazuju  

if is_BIC_used:
    BIC_coeff = ComponentDataBIC.order_oxyfuel_comp1_BIC
else:
    BIC_coeff = None

Util.check_total_fraction(components, title_primer)

#data_nafta , wellstream etc, 250-540 K i 1-190 bara
#oxyfuel itd 260-325 K i 1-100 bara
temperatures = np.arange(230, 540, 1)  
pressures = np.arange(1, 190, 1)      
results = pd.DataFrame(index=pressures, columns=temperatures)
resultsIteration = pd.DataFrame(index=pressures, columns=temperatures)
if toggle_phase_detect:
    excel_rep = pd.DataFrame(columns=["T_K", "P_bar", "FL", "Fv", "Zl", "Zv"])

start_time = time.time() 

if cplusplus:
    temp_list = temperatures.tolist()
    press_list = pressures.tolist()
    cpp_components = Util.convert_to_cpp_components(components)
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
                MT.EOSType.PR,            
                MT.SolveMethod.FSOLVE,
                toggle_phase_detect,
                BIC= BIC_coeff
                ) 
            results.at[Pp, Tt] = Vint = result["V"]
            resultsIteration.at[Pp, Tt] = result["iteration"]
            # if Vint is not None and Vint != -1.0 and Vint != 0:
            #     print(f"V: {Vint}") 
            
            if toggle_phase_detect:
                V_num = result["V"]
                phase = Util.get_phase(V_num)
                if phase == MT.Phase.LIQUID:
                    Vl = 1 
                    Vp = 0
                if phase == MT.Phase.VAPOR:
                    Vl = 0 
                    Vp = 1
                if phase == MT.Phase.VAPORLIQUID: 
                    Vl = 1 - V_num
                    Vp = V_num

                excel_rep.loc[len(excel_rep)] = [Tt-273, Pp, Vl, Vp, result["Zl"], result["Zv"]]

end_time = time.time()  # Kraj mjerenja
elapsed_time = end_time - start_time
print(f"Vrijeme : {elapsed_time:.5f} sek")



#print(results)
#print(resultsIteration)
#results.to_csv("results3.xlsx")
#resultsIteration.to_csv("resultsIT.xlsx")
#V, x, y, method_used, iteration, Zl, Zv
if toggle_phase_detect:
    excel_rep.to_excel("lookup_table_case1_pvt_version.xlsx", index=False, engine="openpyxl")

if not display_iterations:
    results_display = results.astype(float)
else:
    results_display = resultsIteration.astype(int)

Util.display(temperatures, pressures, results_display, title_primer, adjust_display_iterations)




