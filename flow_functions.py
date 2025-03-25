# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 04:30:02 2024

@author: domagoj
"""
from scipy.optimize import fsolve
from numpy.ma.core import log10
import pandas as pd
from scipy.interpolate import griddata
import CoolProp.CoolProp as CP
from numba import jit

def f_Colebrook_White(D, Re, e):
    """
    Solves the Colebrook-White equation using the fsolve function from scipy.optimize.

    Parameters:
        Re (float): Reynolds number
        e (float): Absolute roughness of the pipe, m
        D (float): Diameter of the pipe, m

    Returns:
        f[0] (float): Friction factor
    """
    def f(f):
        return 1 / (f ** 0.5) + 2 * log10(e / (3.7 * D) + 2.51 / (Re * f ** 0.5))
    f0 = 0.01
    f = fsolve(f, f0)
    return f[0]

@jit
def Reynolds(rho, u, D, mu):
    """
    Parameters
    ----------
    rho : float
        density (ρ) of the fluid (SI units: kg/m3).
    u : float
        flow velocity (m/s).
    D : float
        hydraulic diameter (inner diameter of the pipe)
    mu : float
        dynamic viscosity of the fluid (Pa·s or N·s/m2 or kg/(m·s)).

    Returns
    -------
    float
        reynolds number.

    """
    return rho*u*D/mu

@jit
def p_Darcy_Weisbach(v, rho_g, L, f, D):
    """
    Parameters
    ----------
    v : float
        flow velocity (m/s).
    rho_g : float
        density (ρ) of the fluid (SI units: kg/m3).
    L : float
        straight pipe segment length
    f : float
        friction factor.
    D : float
        inner diameter.
    Returns
    -------
    dp : float
        pressure drop, Pa.

    """
    # Inspect changes in viscosity and density, 
    #         as there should be a segment in which the changes are insignificant.
    dp = f*(L/D)*rho_g*(v**2)/2
    return dp

@jit
def dp_table_pure(lookup_table, L, A, d_in, e, p1, T1, qm, nsteps = 10):
  df_dp = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])
  for i in range (nsteps):
    rho_gas = CP.PropsSI('D', 'T', T1, 'P', p1, 'CO2')           # from lookup table
    mu = CP.PropsSI('V', 'T', T1, 'P', p1, 'CO2')                # from lookup table
    qv = qm / rho_gas
    u = qv / A
    Re = Reynolds(rho_gas, u, d_in, mu)
    ff = f_Colebrook_White(d_in, Re, e)
    dP = p_Darcy_Weisbach(v = u, rho_g = rho_gas, L = L/nsteps, f = ff, D = d_in)
    p2 = p1 - dP
    df_dp.loc[i] = [(i+1)*(L/nsteps), p1/1e5, T1-273.15, mu, rho_gas, u, Re, ff, dP/1e5, p2/1e5]
    p1 = p2
  return df_dp

@jit(forceobj=True)
def dp_table(lookup_table, L, A, d_in, e, p1, T1, qm, nsteps = 10):
  points = lookup_table[['p', 't']].values
  df_dp = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])
  for i in range (nsteps):
    rho_gas = griddata(points, lookup_table['rho_g'].values, (p1/1e5, T1-273.15), method='linear')   # from lookup table
    mu = griddata(points, lookup_table['mu_g'].values, (p1/1e5, T1-273.15), method='linear')         # from lookup table
    qv = qm / rho_gas
    u = qv / A
    Re = Reynolds(rho_gas, u, d_in, mu)
    ff = f_Colebrook_White(d_in, Re, e)
    dP = p_Darcy_Weisbach(v = u, rho_g = rho_gas, L = L/nsteps, f = ff, D = d_in)
    p2 = p1 - dP
    df_dp.loc[i] = [(i+1)*(L/nsteps), p1/1e5, T1-273.15, mu, rho_gas, u, Re, ff, dP/1e5, p2/1e5]
    p1 = p2
  return df_dp

