import numpy as np
import math
import pandas as pd
import warnings
import datetime
from numpy.ma.core import log10
from Classes.calc import Calc 
from fastapi.responses import JSONResponse
warnings.filterwarnings("ignore")

def calculate_steps(steps: int, length: int, d_in: float, e: float, p: float, tK: float, qm: float, case: str) -> dict:

    ################################### case handling #################################
    check_case(case)
    df_lookup = None
    is_pure_CO2 = False
    if case.upper()=='CO2':
        is_pure_CO2 = True
    else:
        #TODO: DB
        url = f"https://raw.githubusercontent.com/dvulin/lookup/main/lookup_table_{case}.csv"
        """
        load data from lookup table:
        t (°C), p (bar), FL (liquid fraction), Fv (vapour fraction), rho_L (kg/m3), rho_g (kg/m3),
        mu_L (mPas), mu_g (mPas), Z_oil (compressibility factor), Z_gas (compressibility factor)
        """
        try:
            df_lookup = pd.read_csv(url)
            df_lookup['mu_L']=df_lookup['mu_L']*0.001
            df_lookup['mu_g']=df_lookup['mu_g']*0.001
            print (f'loaded case: {case}')
        except Exception as e:
            print(f"Error loading lookup table: {e}")
            raise e
    ###################################################################################        

    dfi = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])
    TIMEFORMAT = "%H:%M:%S"

    calc_instance = Calc()


    print(f' {p} Pa, {d_in} diameter(m), {steps} steps, {length} length(m), {tK} K, {qm} kg/m3, {e} pipe roughness')
    print('================================================================')           
    print(f'start {datetime.datetime.now().strftime(TIMEFORMAT)}')     


    dfi = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])

    dfi = calc_instance.dp_table_combined(L = length, d_in = d_in, 
                    e = e, p1 = p, T1 = tK, qm = qm, is_pure_CO2 = is_pure_CO2, nsteps = steps, lookup_table = df_lookup)

                        
    dfi['rho_g'] = dfi['rho_g'].astype(float)         # format for plotting
    dfi['mu'] = dfi['mu'].astype(float)               # format for plotting

    print(f'end {datetime.datetime.now().strftime(TIMEFORMAT)}') 
    
    data_dict = dfi.to_dict(orient='records')  # 'records' makes a list of dictionaries
    return data_dict

def check_case(case: str):
    valid_cases = ["case1", "case2", "case3", "CO2"]
    if case not in valid_cases:
        raise ValueError(f"Invalid case: {case}. Must be one of {valid_cases}")
