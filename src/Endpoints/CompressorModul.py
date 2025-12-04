from src.Classes.EOS.ThermodynamicsCalc import CompressorThermo
from src.Classes.UtilClass import Util
from data.testData import ComponentData
import warnings
import datetime
import sys
import os
from typing import List

warnings.filterwarnings("ignore")
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

def compressor_thermodynamics(fractions: List[float], P1: int, P2: int, T1:int, mass_flow: float, isentropic_efficiency: float, polytropic_efficiency: float, full_report: int, polytropic_exponent: float) -> str:
    
    # TIMEFORMAT = "%H:%M:%S"
    # print(f'start {datetime.datetime.now().strftime(TIMEFORMAT)}')  
    
    gass_data = get_and_check_data(fractions)

    ideal_gas_calc = None
    
    if full_report == 1:
        ideal_gas_calc = CompressorThermo.Ideal.calc_ideal_gas_sanity_check(P1, P2, T1, mass_flow, isentropic_efficiency, polytropic_exponent, gass_data)
        #CompressorThermo.Ideal.print_adiabatic_results(P1, T1, P2, mass_flow, isentropic_efficiency, ideal_gas_calc["adiabatic_ideal"])
        #CompressorThermo.Ideal.print_polytropic_results(P1, T1, P2, mass_flow, polytropic_exponent, ideal_gas_calc["polytropic_ideal"])

    real_gas_calc = CompressorThermo.Real.calc_real_gas_thermo(P1, P2, T1, mass_flow, isentropic_efficiency, polytropic_efficiency, gass_data)
    #CompressorThermo.Real.print_adiabatic_results(P1, T1, P2, mass_flow, isentropic_efficiency, ideal_gas_calc["adiabatic_ideal"])
    #CompressorThermo.Real.print_polytropic_results(P1, T1, P2, mass_flow, polytropic_exponent, ideal_gas_calc["polytropic_ideal"])

    # print(f'end {datetime.datetime.now().strftime(TIMEFORMAT)}') 

    match full_report:
        case 1:
            return {
            "ideal": ideal_gas_calc,
            "real": real_gas_calc }
        case 2:
            real_ad    = real_gas_calc["adiabatic_real"]
            real_poly  = real_gas_calc["polytropic_real"] 
            return f"Adiabatic P: {real_ad['P_MW']} MW, W: {real_ad['w_actual']:.3f} kJ/kg | Polytropic P: {real_poly['P_MW']} MW, W: {real_poly['w']:.3f} kJ/kg"
        case _:
            return real_gas_calc


def get_and_check_data(fractions: List[float]):
    gass_data = ComponentData.oxyfuel_comp1

    if len(gass_data) != len(fractions):
        raise ValueError(f"gass_data.count != fractions.count {len(gass_data)}!={len(fractions)}")

    for comp, frac in zip(gass_data, fractions):
        comp.fraction = frac

    Util.check_total_fraction(gass_data, "gass_data")
    return gass_data
                

