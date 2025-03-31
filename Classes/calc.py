import numpy as np
from scipy.optimize import fsolve
from numpy.ma.core import log10
import pandas as pd
from scipy.interpolate import griddata
import CoolProp.CoolProp as CP
from numba import jit

class Calc:
    def __init__(self):
        pass

    def f_Colebrook_White(self, D, Re, e):
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

    def Reynolds(self, rho, u, D, mu):
        """
        Calculates Reynolds number.

        Parameters
        ----------
        rho : float
            density (ρ) of the fluid (SI units: kg/m3).
        u : float
            flow velocity (m/s).
        D : float
            hydraulic diameter (inner diameter of the pipe).
        mu : float
            dynamic viscosity of the fluid (Pa·s or N·s/m2 or kg/(m·s)).

        Returns
        -------
        float
            Reynolds number.
        """
        return rho * u * D / mu


    def p_Darcy_Weisbach(self, v, rho_g, L, f, D):
        """
        Calculates pressure drop using Darcy-Weisbach equation.

        Parameters
        ----------
        v : float
            flow velocity (m/s).
        rho_g : float
            density (ρ) of the fluid (SI units: kg/m3).
        L : float
            straight pipe segment length.
        f : float
            friction factor.
        D : float
            inner diameter.

        Returns
        -------
        dp : float
            pressure drop, Pa.
        """
        dp = f * (L / D) * rho_g * (v ** 2) / 2
        return dp


    def dp_table_combined(self, L, d_in, e, p1, T1, qm, is_pure_CO2, nsteps=10, lookup_table=None):
        
        df_dp = pd.DataFrame(columns=['L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])
        A = 0.25*np.pi*d_in**2

        if not is_pure_CO2:
            points = lookup_table[['p', 't']].values

        for i in range(nsteps):
            if is_pure_CO2:
                # Pure CO2 logic 
                rho_gas = CP.PropsSI('D', 'T', T1, 'P', p1, 'CO2')
                mu = CP.PropsSI('V', 'T', T1, 'P', p1, 'CO2')
            else:
                # Standard lookup table logic
                rho_gas = griddata(points, lookup_table['rho_g'].values, (p1 / 1e5, T1 - 273.15), method='linear')
                mu = griddata(points, lookup_table['mu_g'].values, (p1 / 1e5, T1 - 273.15), method='linear')

            # Common logic for both cases
            qv = qm / rho_gas
            u = qv / A
            Re = self.Reynolds(rho_gas, u, d_in, mu)
            ff = self.f_Colebrook_White(d_in, Re, e)
            dP = self.p_Darcy_Weisbach(v=u, rho_g=rho_gas, L=L / nsteps, f=ff, D=d_in)
            p2 = p1 - dP

            df_dp.loc[i] = [(i+ 1), (i + 1) * (L / nsteps), p1 / 1e5, T1 - 273.15, mu, rho_gas, u, Re, ff, dP / 1e5, p2 / 1e5]
            p1 = p2

        return df_dp

    def Prandtl(self, Cp, mu, lambda_f):
        return Cp * mu / lambda_f

    def dT(self, T0, Ts, Dp, Dt, G, Cp, L, lambda_i, lambda_f, lambda_s, Re_f, Prt, Prf, h):
        """
        Sukhov model for temperature drop in pipeline.
    
        Parameters
        ----------
        T0 : float
            temperature at distance 0 (inlet), K
        Ts : float
            temperature of soil, K
        Dp : float
            outer diameter of pipe, m
        Dt : float
            inner diameter of pipe, m   
        G : float
            mass flow rate, kg/s
        Cp : float
            heat capacity, kJ/kg·K
        L : float
            length of pipeline, km
        lambda_i : float
            thermal conductivity of isolation, W/m·K
        lambda_f : float
            thermal conductivity of fluid, W/m·K
        lambda_s : float
            thermal conductivity of soil, W/m·K
        Re_f : float
            Reynolds number
        Prt : float
            Prandtl number at tube temperature.
        Prf : float
            Prandtl number at fluid temperature.
        h : float
            depth of pipe underground, m

        Returns
        -------
        TL : float
            temperature at distance L (outlet).
        """
        alpha1 = 0.021 * (lambda_f / Dt) * (Re_f ** 0.8) * (Prf ** 0.44) * (Prf / Prt) ** 0.25
        alpha2 = 2 * lambda_s / (0.5 * (Dt + Dp) * np.log(2 * h / (0.5 * (Dt + Dp)) + ((2 * h / (0.5 * (Dt + Dp))) ** 2) - 1))
        K = (1 / (1 / (alpha1 * Dt) + 1 / (alpha2 * Dp) + np.log(Dp / Dt) / (2 * lambda_i))) * 2 / (Dt + Dp)
        TL = T0 + (Ts - T0) * np.exp(-K * np.pi * Dp * L / (G * Cp))
        return TL

    def d_Economic(Q, rho, mu):
        """
            Calculates optimal economic pipe diameter (Zhang, Z.X., Wang, G.X., 
            Massarotto, P., Rudolph, V., 2006. Optimization of pipeline
            transport for CO2 sequestration. 
            Energy Conversion and Management, 47: 702-715.)

        Parameters
        ----------
        Q : float
            mass flow rte kg/s.
        rho : float
            max density of the fluid in a pipe, kg/m3.
        mu : float
            max viscosity of the fluid in a pipe, Pas.
        Returns
        -------
        economic optimum diameter for CO2 transport, float
        """
        return (0.363*((Q/rho)**0.45)*(rho**0.13)*mu**0.025)

    def wall_thickness(P_max, D, S = 483, F = 0.72 , E = 1):
        """
    Parameters
    ----------
    P_max : float
        Max. design pressure (MPa).
    D : float
        Internal diameter (m).
    S : float
        Pipe yield strength (MPa).
    F : float
        Safety (design) factor.
    E : float
        longitudilnal seam joint factor.
    Returns
    -------
    recommended wall thickness, mm
    """
        wt = 0.5*P_max*D/(S*F*E-P_max)
        return wt
