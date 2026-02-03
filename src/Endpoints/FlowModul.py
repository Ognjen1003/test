from src.Classes.Flow.Flow import Flow 
from src.Classes.Flow.LookupTableSingleton import LookupTableSingleton 
from src.EnumsClasses.MethodsAndTypes import CASES
from data.testData import ComponentData
import warnings
import datetime
import sys
import os
import json

warnings.filterwarnings("ignore")
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

#P:pa T:K length:m qm kg/m3 d_in:m e:pipe roughness
def calculate_steps(steps: int, length: int, d_in: float, e: float, p: float, tK: float, qm: float, case: str) -> dict:

    ################################### case handling #################################
    check_case(case)
    
    data = None
    
    if case.upper() == 'CO2':
       case = CASES.CO2
    # elif case.upper() == 'PVT':
    #     case = CASES.PVT
    elif case.upper() == 'OXY1':
            case = CASES.OXY1
            data = ComponentData.oxyfuel_comp1
    elif case.upper() == 'OXY2':
            case = CASES.OXY2
            data = ComponentData.oxyfuel_comp2
    elif case.upper() == 'OXY3':
            case = CASES.OXY3
            data = ComponentData.oxyfuel_comp3
    elif case.upper() == 'CASE1':
        case = CASES.CASE1
        data = LookupTableSingleton.load_lookup_table("case1")
    ###################################################################################        

    TIMEFORMAT = "%H:%M:%S"
    flow_instance = Flow()


    print(f' {p} Pa, {d_in} diameter(m), {steps} steps, {length} length(m), {tK} K, {qm} kg/m3, {e} pipe roughness')
    print('================================================================')           
    print(f'start {datetime.datetime.now().strftime(TIMEFORMAT)}')     

    dfi = flow_instance.dp_table_combined(L = length, d_in = d_in, 
                    e = e, p1 = p, T1 = tK, qm = qm, case=case, nsteps = steps, datasource = data)

                            
    dfi['rho_g'] = dfi['rho_g'].astype(float)         # format for plotting
    dfi['mu'] = dfi['mu'].astype(float)               # format for plotting

    print(f'end {datetime.datetime.now().strftime(TIMEFORMAT)}') 
        
    data_dict = dfi.to_dict(orient='records')  # 'records' makes a list of dictionaries
          

    #zbog ScriptManagera, promijenit u main.py kad tad    
    #return data_dict
    return data_dict

def check_case(case: str):
    valid_cases = ["CASE1", "CO2", "OXY1", "OXY2", "OXY3"]
    if case.upper() not in valid_cases:
        raise ValueError(f"Invalid case: {case}. Must be one of {valid_cases}")


def calculate_steps2(length: float, d_in: float, e: float, p: float, tK: float, qm: float, case: str, fittings=None,
                            step_mode=None, max_step_m=None, virtual_steps=None, viscosity_method=None) -> list[dict]:

    ################################### case handling #################################
    
    check_case(case)

    data = None
    case_upper = case.upper()

    if case_upper == "CO2":
        case_enum = CASES.CO2
    elif case_upper == "OXY1":
        case_enum = CASES.OXY1
        data = ComponentData.oxyfuel_comp1
    elif case_upper == "OXY2":
        case_enum = CASES.OXY2
        data = ComponentData.oxyfuel_comp2
    elif case_upper == "OXY3":
        case_enum = CASES.OXY3
        data = ComponentData.oxyfuel_comp3
    elif case_upper == "CASE1":
        case_enum = CASES.CASE1
        data = LookupTableSingleton.load_lookup_table("case1")
    ################################### case handling #################################

    TIMEFORMAT = "%H:%M:%S"
    flow_instance = Flow()


    print(f"{p} Pa, {d_in} diameter(m), {length} length(m), {tK} K, {qm} kg/s, {e} pipe roughness")
    print("================================================================")
    print(f"start {datetime.datetime.now().strftime(TIMEFORMAT)}")

    # ---- novi poziv (bez nsteps) ----
    dfi_out = flow_instance.dp_table_combined2(
        L=length,
        d_in=d_in,
        e=e,
        p1=p,
        T1=tK,
        qm=qm,
        case=case_enum,
        datasource=data,
        fittings=fittings,
        step_mode=step_mode,
        max_step_m=max_step_m,
        virtual_steps=virtual_steps,
        viscosity_method=viscosity_method,
    )


    dfi = dfi_out[0]
    meta = dfi_out[1] if len(dfi_out) > 1 else None
    # meta samo logiraj
    if meta is not None:
        print(f"dp_table_combined2 meta: {meta}")

    # format for plotting / output
    if "rho_g" in dfi.columns:
        dfi["rho_g"] = dfi["rho_g"].astype(float)
    if "mu" in dfi.columns:
        dfi["mu"] = dfi["mu"].astype(float)

    print(f"end {datetime.datetime.now().strftime(TIMEFORMAT)}")

    data_dict = dfi.to_dict(orient="records")
    return data_dict