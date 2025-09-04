from src.Classes.Flow.Flow import Flow 
from src.Classes.Flow.LookupTableSingleton import LookupTableSingleton 
from src.EnumsClasses.MethodsAndTypes import CASES
import warnings
import datetime
import sys
import os
import json

warnings.filterwarnings("ignore")
sys.path.append(os.path.abspath(os.path.dirname(__file__)))


def calculate_steps(steps: int, length: int, d_in: float, e: float, p: float, tK: float, qm: float, case: str, composition: list) -> dict:

    ################################### case handling #################################
    check_case(case)
    
    df_lookup = None
    
    if case.upper() == 'CO2':
       case = CASES.CO2
    elif case.upper() == 'PVT':
        case = CASES.PVT
        df_lookup = composition
    elif case.upper() == 'CASE1':
        case = CASES.CASE1
        df_lookup = LookupTableSingleton.load_lookup_table("case1")
    ###################################################################################        

    TIMEFORMAT = "%H:%M:%S"
    flow_instance = Flow()


    print(f' {p} Pa, {d_in} diameter(m), {steps} steps, {length} length(m), {tK} K, {qm} kg/m3, {e} pipe roughness')
    print('================================================================')           
    print(f'start {datetime.datetime.now().strftime(TIMEFORMAT)}')     

    dfi = flow_instance.dp_table_combined(L = length, d_in = d_in, 
                    e = e, p1 = p, T1 = tK, qm = qm, case=case, nsteps = steps, lookup_table = df_lookup)

                            
    dfi['rho_g'] = dfi['rho_g'].astype(float)         # format for plotting
    dfi['mu'] = dfi['mu'].astype(float)               # format for plotting

    print(f'end {datetime.datetime.now().strftime(TIMEFORMAT)}') 
        
    data_dict = dfi.to_dict(orient='records')  # 'records' makes a list of dictionaries
          

    #zbog ScriptManagera, promijenit u main.py kad tad    
    #return data_dict
    return json.dumps(data_dict)

def check_case(case: str):
    valid_cases = ["CASE1", "PVT", "CO2"]
    if case.upper() not in valid_cases:
        raise ValueError(f"Invalid case: {case}. Must be one of {valid_cases}")
