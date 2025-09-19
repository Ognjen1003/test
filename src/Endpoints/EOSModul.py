import datetime
from typing import List
from Models.Component import Component
from src.Classes.EOS import RachfordRice, PengRobinsonEOS, SRKEOS, PREnthalpyCalc 
from src.Classes.EOS.RachfordRice import saturation_pressure
from src.EnumsClasses import SolveMethod, EOSType
import math

def perform_eos_calculation(components: List[Component], T: float, P: float, eos_type: EOSType, 
                            method: SolveMethod = SolveMethod.FSOLVE, phase_detect: bool = False, calculate_enthalpy: bool = False) -> dict:
    
    #TIMEFORMAT = "%H:%M:%S"
    #print(f'\n start {datetime.datetime.now().strftime(TIMEFORMAT)}')
    
    eos_models = {
        EOSType.PR: (PengRobinsonEOS, "PengRobinson"),
        EOSType.SRK: (SRKEOS, "SRK")
    }

    if eos_type not in eos_models:
        raise ValueError(f"Unsupported eos_type: {eos_type}")

    eos_class, eos_name = eos_models[eos_type]

    V, x, y, method_used, iteration, Zl, Zv = RachfordRice.solve(eos_class, components, T, P, method=method.value, phase_detect= phase_detect )


    #res = saturation_pressure(eos_class, components, T)

    if phase_detect:
        is_Z_valid = check_Z_validity(Zl, Zv, V)
        if not is_Z_valid:
            raise ValueError(f"Z vrijednosti not valid")


    eos_result = {
        "V": V,
        "x": x,
        "y": y,
        "method": method_used,
        "iteration": iteration,
        "Zl":Zl,
        "Zv":Zv
    }

    # if calculate_enthalpy:
    #     isZValid = check_Z_validity(Zl, Zv, V)
    #     enthalpy_calc = PREnthalpyCalc(components)  
    #     total_enthalpy = enthalpy_calc.get_total_enthalpy(x, y, Zl, Zv, V, T, P)
    #     eos_result["enthalpy"] = total_enthalpy

    #print(f'\n end {datetime.datetime.now().strftime(TIMEFORMAT)}')
    
    return eos_result



def check_Z_validity(Z_l, Z_v, V: float) -> bool:

    def is_valid(Z):
        return Z is not None and isinstance(Z, (int, float)) and Z > 0 and not math.isnan(Z)

    if V == -1:
        return True
    if V == -2 or V == -3 or V == 2:  # samo tekuća faza
        return is_valid(Z_l)
    elif V == -8 or V == -9 or V == 3:  # samo plinska faza
        return is_valid(Z_v)
    else:  # dvije faze
        return is_valid(Z_l) and is_valid(Z_v)


