import datetime
from Classes.Component import Component
from Classes.EOS import RachfordRice, PengRobinsonEOS, SRKEOS
from typing import List



def perform_eos_calculation(components: List[Component], T: float, P: float, z: List[float], eos_type: str) -> dict:

    TIMEFORMAT = "%H:%M:%S"

    print('================================================================')           
    print(f'start {datetime.datetime.now().strftime(TIMEFORMAT)}')     

    """
    components = [
        Component("Methane", 190.6, 4599000, 0.011, 8.07131, 1730.63, 233.426, 19.89, 5.024e-2, 1.269e-5, -11.01e-9),
        Component("Ethane", 305.4, 4872000, 0.099, 8.21201, 1652.57, 229.387, 21.13, 78.60e-3, -11.85e-6, 18.99e-9)
    ]

    T = 300        
    P = 5e6        
    z = [0.5, 0.5] 
    """

    """     V, x, y = RachfordRice.solve(z, PengRobinsonEOS, components, T, P)
    print(f"\n PR EOS — V: {V:.4f}, x: {x}, y: {y}")

    V, x, y = RachfordRice.solve(z, SRKEOS, components, T, P)
    print(f"\n SRK EOS — V: {V:.4f}, x: {x}, y: {y}") 
    """

    if eos_type == "PR":

        V, x, y, method, iteration = RachfordRice.solve(z, PengRobinsonEOS, components, T, P, "fsolve")
        eos_result = {
            "V": V,
            "x": x,
            "y": y,
            "method": method,
            "iteration": iteration,
            "eos_type": "PengRobinson"
        }
        return eos_result

    elif eos_type == "SRK":

        V, x, y, method, iteration = RachfordRice.solve(z, SRKEOS, components, T, P)
        eos_result = {
            "V": V,
            "x": x,
            "y": y,
            "method": method,
            "iteration": iteration,
            "eos_type": "SRK"
        }

        print(f'\n end {datetime.datetime.now().strftime(TIMEFORMAT)}')
        return eos_result
    
    else:
        raise ValueError('{"error": f"Unsupported eos_type: {eos_type}"}')


    

