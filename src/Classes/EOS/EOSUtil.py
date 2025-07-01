from Models.Component import Component
from typing import List
import math

class Calculations:
    """
    Util methods phase equlibrium/Wilson, dew and bubble points - orphans for now
    """
    @staticmethod
    def wilson_K(component: Component, T: float, P: float) -> float:
        Tr = T / component.Tc
        return math.exp(math.log(component.Pc / P) + 5.373 * (1 + component.omega) * (1 - 1 / Tr))
    
    @staticmethod
    def bubble_pressure(components: List[Component], z: List[float], T: float) -> float:
        K = [Calculations.wilson_K(comp, T, 1e5) for comp in components]  # P = 1 bar for  K 
        bubble_P = sum(z[i] * comp.saturation_pressure(T) for i, comp in enumerate(components))
        return bubble_P / sum(z[i] * K[i] for i in range(len(components)))

    @staticmethod
    def dew_pressure(components: List[Component], y: List[float], T: float) -> float:
        K = [Calculations.wilson_K(comp, T, 1e5) for comp in components]  # P = 1 bar for  K 
        dew_P = sum(y[i] * comp.saturation_pressure(T) for i, comp in enumerate(components))
        return dew_P / sum(y[i] * K[i] for i in range(len(components)))