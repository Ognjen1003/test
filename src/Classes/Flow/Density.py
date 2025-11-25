import numpy as np
from typing import Sequence
import src.EnumsClasses.MethodsAndTypes as MT
import src.Classes.UtilClass as Util
from src.Models.Component import Component

class DensityClass:

    @staticmethod 
    # M u kg/mol, P u Pa, T u K, v_shift Peneloux; default: None, ako je v_shift zadan, koristi v_corr = ZRT/P + v_shift.
    def density_from_Z( components: Sequence[Component], T: float, P: float, Z: float, v_shift: float | None = None) -> float:

        if P <= 0 or T <= 0:
            raise ValueError("P i T moraju biti > 0.")
        if Z <= 0:
            raise ValueError("Z mora biti > 0 u fizički smislenom rješenju.")

        x = Util._normalize(np.array([getattr(c, "fraction", 0.0) for c in components], dtype=float))

        # molarna masa smjese (kg/mol)
        M = Util._mw_array_kg_per_mol(components)
        Mbar = float(np.dot(x, M))

        if v_shift is None:
            # standardno iz Z
            return Mbar * P / (float(Z) * MT.CONSTANTS.R * T)
        else:
            # Peneloux korekcija
            v = float(Z) * MT.CONSTANTS.R * T / P          # m^3/mol
            v_corr = v + float(v_shift)
            if v_corr <= 0.0:
                raise ValueError("Korigirani molarni volumen ≤ 0.")
            return Mbar / v_corr
    
    
    @staticmethod
    #V je **molni** udio pare iz flasha.
    def bulk_density(rho_v: float, rho_l: float, V: float) -> float:

        if V < 0.0: V = 0.0
        if V > 1.0: V = 1.0
        if V <= 0.0:
            return float(rho_l)
        if V >= 1.0:
            return float(rho_v)
        if rho_v <= 0 or rho_l <= 0:
            raise ValueError("Gustoce moraju biti > 0")
        return 1.0 / (V / rho_v + (1.0 - V) / rho_l)
