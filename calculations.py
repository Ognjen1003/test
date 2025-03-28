import numpy as np
import math
import pandas as pd
import warnings
import os
import datetime
#from flow_functions import dp_table
from numpy.ma.core import log10
from Classes.calc import Calc 
warnings.filterwarnings("ignore")

def calculate_steps(steps: int, length: int, d_in: float, e: float, p: float, tK: float, qm: float, case: str) -> float:
    
    ################################### case handling #################################
    url = f"https://raw.githubusercontent.com/dvulin/lookup/main/lookup_table_{case}.csv"

    """
    load data from lookup table:
    t (°C), p (bar), FL (liquid fraction), Fv (vapour fraction), rho_L (kg/m3), rho_g (kg/m3),
    mu_L (mPas), mu_g (mPas), Z_oil (compressibility factor), Z_gas (compressibility factor)
    """
    df_lookup = pd.read_csv(url)
    df_lookup['mu_L']=df_lookup['mu_L']*0.001
    df_lookup['mu_g']=df_lookup['mu_g']*0.001
    print (f'loaded case: {case}')
    ###################################################################################        

    dfi = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])
    TIMEFORMAT = "%H:%M:%S"
    p_out = 0.0                         #moja varijabla
    area = 0.25*np.pi*d_in**2
    
    Dexp = 0.315925 

    parameter_sensitivity = {}
    d_in_sens = {}
    e_sens = {}      

    calc_instance = Calc()


    print(f'A {area}, {p} Pa, {d_in} diameter(m), {steps} steps, {length} length(m), {tK} K, {qm} kg/m3, {e} pipe roughness')
    print('================================================================')           
    print(f'start {datetime.datetime.now().strftime(TIMEFORMAT)}')     


    dfi = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])
    dfi = calc_instance.dp_table(lookup_table = df_lookup, L = length, A = area, d_in = d_in, e = e, p1 = p, T1=tK, qm = qm, nsteps = steps)
                        
    dfi['rho_g'] = dfi['rho_g'].astype(float)         # format for plotting
    dfi['mu'] = dfi['mu'].astype(float)               # format for plotting
               
    e_sens[e] = dfi.copy()
    d_in_sens[d_in] = e_sens.copy()
    parameter_sensitivity[p] = d_in_sens.copy()

    p2 = parameter_sensitivity[p][d_in][e]['p2'].iloc[-1]

    try:
        p_out = p2
        if math.isnan(float(p2)):
            raise ValueError("The calculated value is NaN")
        print(f"================================ pad tlaka: {p_out}")
    except ValueError as e:
        print(f"Error: {e}")
        raise ValueError(e)


    print(f'end {datetime.datetime.now().strftime(TIMEFORMAT)}') 
    
    return p_out