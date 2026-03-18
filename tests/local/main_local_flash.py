from src.Endpoints.EOSModul import perform_eos_calculation
from src.Classes.UtilClass import Util
from data.testData import ComponentData
from data.testDataBIC import ComponentDataBIC
import src.EnumsClasses.MethodsAndTypes as MT
import pandas as pd
import numpy as np
import time, sys, os
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

bin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
sys.path.insert(0, bin_path)

print("Current dir:", os.getcwd())
print("sys.path:", sys.path)

import eos_cpp


# u bin imas pyd file, to je dll sa zvanje iz cpp
cplusplus = False
display_iterations = False
adjust_display_iterations = False
toggle_phase_detect = False
is_BIC_used = False


title_primer = "oxyfuel_comp1"
components = ComponentData.data_components

if is_BIC_used:
    BIC_coeff = ComponentDataBIC.order_oxyfuel_comp1_BIC
else:
    BIC_coeff = None


def solve_one_temperature(args):
    Tt, pressures, components, toggle_phase_detect, BIC_coeff = args

    row_results = []
    excel_rows = []

    for Pp in pressures:
        result = perform_eos_calculation(
            components,
            Tt,
            Pp,
            MT.EOSType.SRK,
            MT.SolveMethod.FSOLVE,
            toggle_phase_detect,
            BIC=BIC_coeff
        )

        row_results.append((Tt, Pp, result["V"], result["iteration"]))

        if toggle_phase_detect:
            V_num = result["V"]
            phase = Util.get_phase(V_num)

            if phase == MT.Phase.LIQUID:
                Vl = 1
                Vp = 0
            elif phase == MT.Phase.VAPOR:
                Vl = 0
                Vp = 1
            elif phase == MT.Phase.VAPORLIQUID:
                Vl = 1 - V_num
                Vp = V_num
            else:
                Vl = None
                Vp = None

            excel_rows.append((Tt - 273, Pp, Vl, Vp, result["Zl"], result["Zv"]))

    return row_results, excel_rows


def main():
    Util.check_total_fraction(components, title_primer)

    # data_nafta , wellstream etc, 250-540 K i 1-190 bara
    temperatures = np.arange(230, 540, 1)
    pressures = np.arange(1, 190, 1)

    results = pd.DataFrame(index=pressures, columns=temperatures)
    resultsIteration = pd.DataFrame(index=pressures, columns=temperatures)

    if toggle_phase_detect:
        excel_rep = pd.DataFrame(columns=["T_K", "P_bar", "FL", "Fv", "Zl", "Zv"])
    else:
        excel_rep = None

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
        task_args = [
            (Tt, pressures, components, toggle_phase_detect, BIC_coeff)
            for Tt in temperatures
        ]

        max_workers = max(1, os.cpu_count() - 1)

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            for row_results, excel_rows in executor.map(solve_one_temperature, task_args):
                for Tt, Pp, Vval, iteration in row_results:
                    results.at[Pp, Tt] = Vval
                    resultsIteration.at[Pp, Tt] = iteration

                if toggle_phase_detect and excel_rows:
                    for row in excel_rows:
                        excel_rep.loc[len(excel_rep)] = row

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Vrijeme : {elapsed_time:.5f} sek")

    if toggle_phase_detect:
        excel_rep.to_excel("lookup_table_case1_pvt_version.xlsx", index=False, engine="openpyxl")

    if not display_iterations:
        results_display = results.astype(float)
    else:
        results_display = resultsIteration.astype(int)

    Util.display(temperatures, pressures, results_display, title_primer, adjust_display_iterations)


if __name__ == "__main__":
    mp.freeze_support()
    main()