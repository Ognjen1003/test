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
import sys, os, math




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
    def dp_table_combined2(self, L, d_in, e, p1, T1, qm, case, datasource=None, fittings=None, step_mode=MT.StepMode.BY_FITTINGS, 
                                max_step_m=5000.0, virtual_steps=None, viscosity_method=None):
        
        breakpoints, K_at, d_in_new_at, K_ref_at, meta = self.fitting_validator(
            L=L,
            fittings=fittings,
            step_mode=step_mode,
            max_step_m=max_step_m,
            virtual_steps=virtual_steps
        )

        df_dp = pd.DataFrame(columns=[
            "step", "L", "p1", "t", "mu", "rho_g", "u", "Re", "ff",
            "K_step", "dp_fric", "dp_minor", "dp_total", "p2",
            "d_in"  # korisno za debug kad ima promjena promjera
        ])

        p1 = float(p1)  # Pa
        d_in = float(d_in)  # trenutni promjer (m)

        for i in range(len(breakpoints) - 1):
            x0 = breakpoints[i]
            x1 = breakpoints[i + 1]
            L_seg = x1 - x0
            if L_seg <= 0:
                continue

            A = 0.25 * np.pi * d_in**2

            if case == MT.CASES.CO2:
                rho = CP.PropsSI('D', 'T', T1, 'P', p1, 'CO2')
                mu  = CP.PropsSI('V', 'T', T1, 'P', p1, 'CO2')

            elif case in (MT.CASES.OXY1, MT.CASES.OXY2, MT.CASES.OXY3):
                if datasource is None:
                    raise ValueError("Composition can't be empty for OXY alike cases")
                composition = datasource
                p_bar = p1 / 100000.0
                rho = self.get_rho(p1, T1, composition, p_bar)
                mu  = self.get_viscosity(composition, T1, rho)

            elif case == MT.CASES.CASE1:
                points = datasource[['p', 't']].values
                rho = griddata(points, datasource['rho_g'].values, (p1 / 1e5, T1 - 273.15), method='linear')
                mu  = griddata(points, datasource['mu_g'].values,  (p1 / 1e5, T1 - 273.15), method='linear')

            qv = qm / rho
            u = qv / A
            Re = self.Reynolds(rho, u, d_in, mu)
            ff = self.f_Colebrook_White(d_in, Re, e)

            dP_fric = self.p_Darcy_Weisbach(v=u, rho_g=rho, L=L_seg, f=ff, D=d_in)
            p2 = p1 - dP_fric

            at_key = round(float(x1), 6)
            K_step = float(K_at.get(at_key, 0.0))
            d_in_next = d_in_new_at.get(at_key, None)
            dP_minor = 0.0

            if K_step != 0.0:
                K_ref = (K_ref_at.get(at_key, None) or "UPSTREAM").upper()

                # odabir referentnog promjera za brzinu u K * (rho*u^2/2)
                # - UPSTREAM: koristi trenutni d_in (prije promjene)
                # - DOWNSTREAM: koristi d_in_next (nakon promjene)
                # - MIN_DIAMETER: koristi min(d_in, d_in_next)
                if d_in_next is not None and K_ref in ("DOWNSTREAM", "MIN_DIAMETER"):
                    if K_ref == "DOWNSTREAM":
                        d_ref = float(d_in_next)
                    else:  # MIN_DIAMETER
                        d_ref = float(min(d_in, d_in_next))
                else:
                    d_ref = float(d_in)

                A_ref = 0.25 * np.pi * d_ref**2
                u_ref = qv / A_ref

                dP_minor = K_step * (0.5 * rho * u_ref**2)   # Pa
                p2 = p2 - dP_minor

            dP_total = dP_fric + dP_minor
            df_dp.loc[i] = [
                i + 1,
                x1,                 # kumulativna duljina (m)
                p1 / 1e5,            # bar (kao prije)
                T1 - 273.15,
                mu,
                rho,
                u,
                Re,
                ff,
                K_step,
                dP_fric / 1e5,       # bar
                dP_minor / 1e5,      # bar
                dP_total / 1e5,      # bar
                p2 / 1e5,            # bar
                d_in
            ]

            if at_key in d_in_new_at:
                d_in = float(d_in_new_at[at_key])

            # idući korak
            p1 = p2

        return df_dp, meta

    
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
            
    def get_viscosity(self, composition, T1, rho, viscosity_method = None):
        if viscosity_method is None or viscosity_method is MT.ViscosityMethod.AUTO or viscosity_method is MT.ViscosityMethod.WILKE:
            viscosity = V.ViscosityClass.gas_viscosity_wilke(composition, T1)
        else: 
            raise ValueError("Not implemented") 
        #viscosity = V.ViscosityClass.viscosity_lbc(composition, T1, rho)
        return viscosity

    import math

    def fitting_validator(self, *, L: float, fittings: list | None, step_mode, max_step_m: float | None = None, virtual_steps: int | None = None):
        if L is None or float(L) <= 0:
            raise ValueError("L mora biti > 0")
        L = float(L)

        if fittings is None:
            fittings = []

        # ---------------------------------------------------
        # 1) Virtual steps: ravna cijev, N segmenata
        # ---------------------------------------------------
        if virtual_steps is not None:
            if not isinstance(virtual_steps, int) or virtual_steps <= 0:
                raise ValueError("virtual_steps mora biti pozitivan int")

            step_len = L / virtual_steps
            breakpoints = [i * step_len for i in range(0, virtual_steps + 1)]
            breakpoints[-1] = L

            return breakpoints, {}, {}, {}, {
                "mode": "VIRTUAL_STEPS",
                "virtual_steps": virtual_steps,
                "fittings_count": 0
            }

        # ---------------------------------------------------
        # 2) Clean fittings: dict sa at_m, K, d_in_new, K_ref
        # ---------------------------------------------------
        clean_fittings = []

        for f in fittings:
            at_m = getattr(f, "at_m", None) if not isinstance(f, dict) else f.get("at_m")
            Kval = getattr(f, "K", None) if not isinstance(f, dict) else f.get("K")
            d_in_new = getattr(f, "d_in_new", None) if not isinstance(f, dict) else f.get("d_in_new")
            K_ref = getattr(f, "K_ref", None) if not isinstance(f, dict) else f.get("K_ref")

            if at_m is None:
                raise ValueError("Fitting mora imati 'at_m'")

            if Kval is None and d_in_new is None:
                raise ValueError(f"Fitting na at_m={at_m} mora imati barem 'K' ili 'd_in_new'")

            at_m = float(at_m)
            if at_m < 0 or at_m > L:
                raise ValueError(f"Fitting at_m={at_m} je izvan raspona [0, L={L}]")

            if Kval is not None:
                Kval = float(Kval)

            if d_in_new is not None:
                d_in_new = float(d_in_new)
                if d_in_new <= 0:
                    raise ValueError(f"d_in_new mora biti > 0 (at_m={at_m})")

            if K_ref is None:
                K_ref = "UPSTREAM"
            else:
                K_ref = str(K_ref).upper()
                if K_ref not in ("UPSTREAM", "DOWNSTREAM", "MIN_DIAMETER"):
                    raise ValueError("K_ref mora biti: UPSTREAM, DOWNSTREAM ili MIN_DIAMETER")

            clean_fittings.append({
                "at_m": at_m,
                "K": Kval,
                "d_in_new": d_in_new,
                "K_ref": K_ref
            })

        clean_fittings.sort(key=lambda x: x["at_m"])

        # ---------------------------------------------------
        # 3) Mapiranje na breakpoint ključ (round)
        # ---------------------------------------------------
        K_at = {}
        d_in_new_at = {}
        K_ref_at = {}

        for f in clean_fittings:
            at_key = round(f["at_m"], 6)

            if f["K"] is not None:
                K_at[at_key] = K_at.get(at_key, 0.0) + f["K"]

            if f["d_in_new"] is not None:
                d_in_new_at[at_key] = f["d_in_new"]

            # ako već postoji, ostavi prvi (ili zadnji – po želji)
            if at_key not in K_ref_at:
                K_ref_at[at_key] = f["K_ref"]

        # breakpointovi koje moramo imati (spajamo i K i promjere)
        fitting_points = sorted(set(list(K_at.keys()) + list(d_in_new_at.keys())))

        # ---------------------------------------------------
        # 4) Generiranje breakpoints
        # ---------------------------------------------------
        # Default: BY_FITTINGS, a max_step_m je opcionalni "refine"
        if getattr(step_mode, "name", str(step_mode)) == "BY_FITTINGS":
            breakpoints = [0.0] + fitting_points + [L]

            if max_step_m is not None and float(max_step_m) > 0:
                max_step_m = float(max_step_m)
                refined = [0.0]
                for a, b in zip(breakpoints[:-1], breakpoints[1:]):
                    seg_len = b - a
                    if seg_len <= max_step_m:
                        refined.append(b)
                    else:
                        n = int(math.ceil(seg_len / max_step_m))
                        dl = seg_len / n
                        for k in range(1, n + 1):
                            refined.append(a + k * dl)
                refined[-1] = L
                breakpoints = refined

        else:
            # FIXED: korak je max_step_m
            if max_step_m is None or float(max_step_m) <= 0:
                raise ValueError("Za FIXED mode treba max_step_m > 0")

            max_step_m = float(max_step_m)
            n = int(math.ceil(L / max_step_m))
            dl = L / n
            breakpoints = [i * dl for i in range(n + 1)]
            breakpoints[-1] = L

        # ukloni duplikate
        bp = []
        for x in breakpoints:
            x = float(x)
            if not bp or abs(x - bp[-1]) > 1e-9:
                bp.append(x)
        breakpoints = bp

        if breakpoints[0] != 0.0 or abs(breakpoints[-1] - L) > 1e-6:
            raise ValueError("Breakpointovi nisu dobro generirani (ne počinju s 0 ili ne završavaju s L).")

        return breakpoints, K_at, d_in_new_at, K_ref_at, {
            "mode": step_mode.value if hasattr(step_mode, "value") else str(step_mode),
            "virtual_steps": None,
            "fittings_count": len(clean_fittings),
            "breakpoints_count": len(breakpoints)
        }

