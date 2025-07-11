import pandas as pd
import warnings
import datetime
from src.Classes.Flow.Flow import Flow 
from src.Classes.Flow.LookupTableSingleton import LookupTableSingleton 
import sys
import os
import json

warnings.filterwarnings("ignore")
sys.path.append(os.path.abspath(os.path.dirname(__file__)))


def calculate_steps(steps: int, length: int, d_in: float, e: float, p: float, tK: float, qm: float, case: str) -> dict:

    ################################### case handling #################################
    check_case(case)
    df_lookup = None
    is_pure_CO2 = False
    if case.upper()=='CO2':
        is_pure_CO2 = True
    else:
        #TODO: DB
        df_lookup = LookupTableSingleton.load_lookup_table("case1")
    ###################################################################################        

    dfi = pd.DataFrame(columns=['step', 'L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])
    TIMEFORMAT = "%H:%M:%S"

    flow_instance = Flow()


    print(f' {p} Pa, {d_in} diameter(m), {steps} steps, {length} length(m), {tK} K, {qm} kg/m3, {e} pipe roughness')
    print('================================================================')           
    print(f'start {datetime.datetime.now().strftime(TIMEFORMAT)}')     


    dfi = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])

    dfi = flow_instance.dp_table_combined(L = length, d_in = d_in, 
                    e = e, p1 = p, T1 = tK, qm = qm, is_pure_CO2 = is_pure_CO2, nsteps = steps, lookup_table = df_lookup)

                        
    dfi['rho_g'] = dfi['rho_g'].astype(float)         # format for plotting
    dfi['mu'] = dfi['mu'].astype(float)               # format for plotting

    print(f'end {datetime.datetime.now().strftime(TIMEFORMAT)}') 
    
    data_dict = dfi.to_dict(orient='records')  # 'records' makes a list of dictionaries

    #zbog ScriptManagera, promijenit u main.py kad tad    
    #return data_dict
    return json.dumps(data_dict)

def check_case(case: str):
    valid_cases = ["case1", "case2", "case3", "CO2"]
    if case not in valid_cases:
        raise ValueError(f"Invalid case: {case}. Must be one of {valid_cases}")
