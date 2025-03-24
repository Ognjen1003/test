from numpy.ma.core import log10
import numpy as np
import pandas as pd
import warnings
import os
warnings.filterwarnings("ignore")


def calculate_division(steps: int, length: int) -> float:
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
    L = 40*1000                           #@param L (float) = 40*1000       # pipe length, m
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
    nsteps = 40                         # @param nsteps (int) = 25  # number of steps
    boosters = False                     # @param {type:"boolean"}


    p_l, t_l, dp_l, rho_g_l = [], [], [], []
    dfi = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])

    d_in = np.array([241.1])/1000             # m
    wthick = np.array([16, 16 ,16, 16, 16])                             #
    p_in = (np.array([40, 50])*1e5)                         # Pa
    e_i = np.array([0.05, 0.025])/1000                # m
    T = np.array([40])+273.15                       # K
    Q = np.array([752])*1e6/(365*24*3600)   

    import datetime
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
                now = datetime.datetime.now()
                print (f'+- p_in = {p1} Pa')
                for D in d_in:
                    now = datetime.datetime.now()
                    print (f'   +--- {now.strftime(timeformat)}:  d_in = {D} m')
                    dfi = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])
                    e_sens = {}
                    for e in e_i:
                        p_1 = p1
                        now = datetime.datetime.now()
                        print (f'         +--- {now.strftime(timeformat)}:  roughness = {e} m')
                        dfi = dp_table(lookup_table = df_lookup, 
                                    L = L, A = A, d_in = D, 
                                    e = e, p1 = p1, T1=T1, 
                                    qm = qm, nsteps = nsteps)
                        dfi['rho_g'] = dfi['rho_g'].astype(float)         # format for plotting
                        dfi['mu'] = dfi['mu'].astype(float)               # format for plotting
                        e_sens[e] = dfi.copy()
                    d_in_sens[D] = e_sens.copy()
                parameter_sensitivity[p1] = d_in_sens.copy()
            
            df_summary = pd.DataFrame(columns=['p_in', 'd_in', 'e', 'p_out'])
            for p1 in p_in: 
                for D in d_in:
                    for e in e_i: 
                        p2 = parameter_sensitivity[p1][D][e]['p2'].iloc[-1]
                        print(f'p2: {p2}')
                        new_row = {'p_in': p1/1e5, 'd_in':D, 'e': e, 'p_out': p2}
                        df_summary = df_summary._append(new_row, ignore_index = True)
            directory = f'case_{int(qm/(1e6 / (365*24*3600)))}/'
            if not os.path.exists(directory):
                os.makedirs(directory)
            df_summary.to_csv(f'{directory}{case}_summary_nsteps_{nsteps}_t_{int(T1-273.15)}.txt')
            print (f'saved to: {directory}{case}_summary_nsteps_{nsteps}_t_{int(T1-273.15)}.txt')
        
            with open(f'{directory}{case}_sensitivity_t_{int(T1-273.15)}.json', 'w') as f:
                for p_i, d_in_dict in parameter_sensitivity.items():
                    for din, e_sens_dict in d_in_dict.items():
                        for ei, df in e_sens_dict.items():
                            df_json = df.to_json(orient='split')
                            f.write(f'{p_i},{din},{ei},{df_json}\n')

        end_t = datetime.datetime.now()
        runtime = end_t - start_t
        print (f'ended at {end_t}, runtime: {runtime.total_seconds():.1f} seconds' )
    

    return 0.0