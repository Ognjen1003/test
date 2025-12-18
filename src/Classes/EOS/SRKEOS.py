import math
from .EOS import EOSBase
import src.Models.Component as Component
from typing import List

class SRKEOS(EOSBase):
    """
    Implementacija SRK jednadžbe stanja.
    """
    def alpha(self, comp: Component) -> float:
        Tr = self.T / comp.Tc
        m = 0.48 + 1.574 * comp.omega - 0.176 * comp.omega**2
        return (1 + m * (1 - math.sqrt(Tr)))**2

    def a_formula(self, comp: Component) -> float:
        return 0.42748 * self.R**2 * comp.Tc**2 / comp.Pc

    def b_formula(self, comp: Component) -> float:
        return 0.08664 * self.R * comp.Tc / comp.Pc

    def get_coefficients(self, A: float, B: float) -> List[float]:
        return [1, -1, A - B - B**2, -A * B]