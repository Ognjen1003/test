from src.Endpoints.EOSModul import perform_eos_calculation
import src.EnumsClasses.MethodsAndTypes as MT
import src.Classes.Flow.Density as D
import src.Classes.Flow.Viscosity as V
import src.Classes.UtilClass as Util
import pandas as pd
from numpy.ma.core import log10
import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import griddata
import CoolProp.CoolProp as CP
from numba import jit
import sys, os




bin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\\..', '..', 'bin'))
sys.path.insert(0, bin_path)  # stavim ga kao prvi prioritet

import eos_cpp 

class Flow:
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

    #@jit(forceobj=True)
    def dp_table_combined(self, L, d_in, e, p1, T1, qm, case, nsteps=10, datasource=None):
        
        df_dp = pd.DataFrame(columns=['step', 'L', 'p1', 't', 'mu', 'rho_g', 'u', 'Re', 'ff', 'dp', 'p2'])
        A = 0.25*np.pi*d_in**2

        if case == MT.CASES.CASE1:
            points = datasource[['p', 't']].values

        for i in range(nsteps):
            if case == MT.CASES.CO2:
                # Pure CO2 logic (from dp_table_pure)
                rho = CP.PropsSI('D', 'T', T1, 'P', p1, 'CO2')
                mu = CP.PropsSI('V', 'T', T1, 'P', p1, 'CO2')
            elif case == MT.CASES.OXY1 or case == MT.CASES.OXY2 or case == MT.CASES.OXY3:
                composition = datasource
                if datasource is None:
                    raise ValueError("Composition can't be empty for OXY alike cases")
                p_bar = p1/100000
                rho = self.get_rho(p1, T1, composition, p_bar)
                mu = self.get_viscosity(composition, T1, rho)
                
            elif case == MT.CASES.CASE1:
                # Standard lookup table logic (from dp_table)
                rho = griddata(points, datasource['rho_g'].values, (p1 / 1e5, T1 - 273.15), method='linear')
                mu = griddata(points, datasource['mu_g'].values, (p1 / 1e5, T1 - 273.15), method='linear')

            # Common logic for both cases
            qv = qm / rho
            u = qv / A
            Re = self.Reynolds(rho, u, d_in, mu)
            ff = self.f_Colebrook_White(d_in, Re, e)
            dP = self.p_Darcy_Weisbach(v=u, rho_g=rho, L=L / nsteps, f=ff, D=d_in)
            p2 = p1 - dP

            df_dp.loc[i] = [(i+ 1), (i + 1) * (L / nsteps), p1 / 1e5, T1 - 273.15, mu, rho, u, Re, ff, dP / 1e5, p2 / 1e5]
            p1 = p2

        return df_dp

    def get_rho(self, p1, T1, composition, p_bar):

        result = perform_eos_calculation(composition, T1, p_bar, MT.EOSType.PR, MT.SolveMethod.FSOLVE, True) 
                
        V = result["V"]
        Zv = result["Zv"]
        Zl = result["Zl"]

        actual_phase = Util.Util.get_phase(V)

        if actual_phase == MT.Phase.LIQUID:                     
            rho = D.DensityClass.density_from_Z(composition, T1, p1, Zl)
        elif actual_phase == MT.Phase.VAPOR:                     
            rho = D.DensityClass.density_from_Z(composition, T1, p1, Zv)
        elif actual_phase == MT.Phase.VAPORLIQUID:
            rho_liquid = D.DensityClass.density_from_Z(composition, T1, p1, Zl)
            rho_vapor = D.DensityClass.density_from_Z(composition, T1, p1, Zv)
            rho = D.DensityClass.bulk_density(rho_vapor, rho_liquid, V)
        else:
            raise ValueError("phase problem")
        
        return rho
            
    def get_viscosity(self, composition, T1, rho):
        viscosity = V.ViscosityClass.gas_viscosity_wilke(composition, T1)
        #viscosity = V.ViscosityClass.viscosity_lbc(composition, T1, rho)
        return viscosity