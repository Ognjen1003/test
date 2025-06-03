import datetime
from typing import List
from src.Classes.Component import Component
from src.Classes.EOS import RachfordRice, PengRobinsonEOS, SRKEOS
from src.EnumsClasses import SolveMethod, EOSType

def perform_eos_calculation(components: List[Component], T: float, P: float, eos_type: EOSType, method: SolveMethod = SolveMethod.FSOLVE) -> dict:
    
    #TIMEFORMAT = "%H:%M:%S"
    #print(f'\n start {datetime.datetime.now().strftime(TIMEFORMAT)}')
    
    eos_models = {
        EOSType.PR: (PengRobinsonEOS, "PengRobinson"),
        EOSType.SRK: (SRKEOS, "SRK")
    }

    if eos_type not in eos_models:
        raise ValueError(f"Unsupported eos_type: {eos_type}")

    eos_class, eos_name = eos_models[eos_type]

    V, x, y, method_used, iteration = RachfordRice.solve(eos_class, components, T, P, method=method.value)

    eos_result = {
        "V": V,
        "x": x,
        "y": y,
        "method": method_used,
        "iteration": iteration,
        "eos_type": eos_name
    }

    #print(f'\n end {datetime.datetime.now().strftime(TIMEFORMAT)}')
    
    return eos_result

