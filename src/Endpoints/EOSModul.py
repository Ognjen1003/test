import src.EnumsClasses.MethodsAndTypes as MT
from src.Models.Component import Component
from src.Classes.EOS import RachfordRice, PengRobinsonEOS, SRKEOS 
from src.EnumsClasses import SolveMethod, EOSType
from src.Classes.UtilClass import Util
import datetime 
import numpy as np
from typing import List
import math

def perform_eos_calculation(components: List[Component], T_K: float, P_bar: float, eos_type: EOSType, 
                            method: SolveMethod = SolveMethod.FSOLVE, phase_detect: bool = False, BIC: np.ndarray | None = None ) -> dict:
    
    #TIMEFORMAT = "%H:%M:%S"
    #print(f'\n start {datetime.datetime.now().strftime(TIMEFORMAT)}')

    eos_models = {
        EOSType.PR: (PengRobinsonEOS, "PengRobinson"),
        EOSType.SRK: (SRKEOS, "SRK")
    }

    if eos_type not in eos_models:
        raise ValueError(f"Unsupported eos_type: {eos_type}")

    eos_class, eos_name = eos_models[eos_type]

    V, x, y, method_used, iteration, Zl, Zv = RachfordRice.solve(eos_class, components, T_K, P_bar, method=method.value, phase_detect= phase_detect, k_ij= BIC )


    #res = saturation_pressure(eos_class, components, T)

    if phase_detect:
        is_Z_valid = check_Z_validity(Zl, Zv, V)
        if not is_Z_valid:
            raise ValueError(f"Z vrijednosti not valid P_bar:{P_bar} - T:{T_K}")


    eos_result = {
        "V": V,
        "x": x,
        "y": y,
        "method": method_used,
        "iteration": iteration,
        "Zl":Zl,
        "Zv":Zv
    }


    #print(f'\n end {datetime.datetime.now().strftime(TIMEFORMAT)}')
    
    return eos_result



def check_Z_validity(Z_l, Z_v, V: float) -> bool:

    def is_valid(Z):
        return Z is not None and isinstance(Z, (int, float)) and Z > 0 and not math.isnan(Z)

    phase = Util.get_phase(V)

    if V == -1:
        return True
    if phase == MT.Phase.LIQUID:
        return is_valid(Z_l)
    elif phase == MT.Phase.VAPOR:  # samo plinska faza
        return is_valid(Z_v)
    else:  # dvije faze
        return is_valid(Z_l) and is_valid(Z_v)


