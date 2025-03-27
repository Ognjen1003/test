import numpy as np
import pandas as pd
import warnings
import os
import datetime
from flow_functions import dp_table
from other_functions import *
from numpy.ma.core import log10
warnings.filterwarnings("ignore")

def calculate_division(steps: int, length: int, d_in: float, roughness: float, p: float, T: float, qm: float) -> float:
    

    pure_CO2 = False                     # @param {type:"boolean"}

    case = 'case1'
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


    p1 =50*1e5                            #@param p1 (float) = 40*1e5       # p, Pa
    T1 = 273.15+50                        #@param T1 (float) = 273.15+60    # T, K
    L = length                           #@param L (float) = 40*1000       # pipe length, m
    e = 0.0457 / 1000                     #@param e (float) = 0.0457 / 1000 # pipe roughness, m
    D = 0.315925                          #@param D (float) = 0.325   # m

    """
    za odabir posluzilo:
    Peletiri 2018, Table 1 - 0.325 m?
    u literaturi za CO2 iz postrojenja za etanol, rasponi su 4 do 8 inch za priključne cjevovode
    """
    A = 0.25*np.pi*D**2
    qm = 650                             # @param qm (float) = 700   # ktpa
    qm = qm * 1e6 / (365*24*3600)        # mass flow rate (kg/s) of CO2 stream = 700 kt/y = 22.19685438863521 kg/s
    nsteps = steps                         # @param nsteps (int) = 25  # number of steps
    boosters = False                     # @param {type:"boolean"}


    p_l, t_l, dp_l, rho_g_l = [], [], [], []
    dfi = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])

    d_in = np.array([241.1])/1000             # m
    wthick = np.array([16, 16 ,16, 16, 16])                             #
    p_in = (np.array([40])*1e5)                         # Pa
    e_i = np.array([0.05])/1000                # m
    T = np.array([40])+273.15                       # K
    Q = np.array([750])*1e6/(365*24*3600)
    p_last = 0.0   #moja varijabla


    timeformat = "%H:%M:%S"
    start_t = datetime.datetime.now()
    parameter_sensitivity = {}
    for qm in Q:
        print ('----------------------------------------------------------------------')
        print (f'qm = {int(qm/(1e6 / (365*24*3600)))} ktpa')
        for T1 in T:
            print (f'T = {T1-273.15}°C')
            for p1 in p_in:
                d_in_sens = {}
                print (f'+- p_in = {p1} Pa')
                for D in d_in:
                    print (f'd_in = {D} m')
                    dfi = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])
                    e_sens = {}
                    for e in e_i:
                        print (f'+ begin --- {datetime.datetime.now().strftime(timeformat)}: roughness = {e} m')
                        dfi = dp_table(lookup_table = df_lookup, 
                                    L = L, A = A, d_in = D, 
                                    e = e, p1 = p1, T1=T1, 
                                    qm = qm, nsteps = nsteps)
                        print (f'+ end --- {datetime.datetime.now().strftime(timeformat)}: roughness = {e} m')
                        dfi['rho_g'] = dfi['rho_g'].astype(float)         # format for plotting
                        dfi['mu'] = dfi['mu'].astype(float)               # format for plotting
                        e_sens[e] = dfi.copy()
                    d_in_sens[D] = e_sens.copy()
                parameter_sensitivity[p1] = d_in_sens.copy()
            

            for p1 in p_in: 
                for D in d_in:
                    for e in e_i: 
                        p2 = parameter_sensitivity[p1][D][e]['p2'].iloc[-1]
                        p_last = p2

    print(f"================================ pad tlaka: {p_last}")
    return p_last