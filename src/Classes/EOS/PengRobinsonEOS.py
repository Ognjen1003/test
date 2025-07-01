import math
from .EOS import EOSBase
import Models.Component as Component
from typing import List


class PengRobinsonEOS(EOSBase):
    """
    Implementacija PR jednadžbe stanja.
    """
    
    def alpha(self, comp: Component) -> float:
        Tr = self.T / comp.Tc
        m = 0.37464 + 1.54226 * comp.omega - 0.26992 * comp.omega**2
        return (1 + m * (1 - math.sqrt(Tr)))**2

    def a_formula(self, comp: Component) -> float:
        return 0.45724 * self.R**2 * comp.Tc**2 / comp.Pc

    def b_formula(self, comp: Component) -> float:
        return 0.07780 * self.R * comp.Tc / comp.Pc

    def get_coefficients(self, A: float, B: float) -> List[float]:
        return [1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]


