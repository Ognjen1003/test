import math
import numpy as np
from .EOS import EOSBase
import src.Models.Component as Component
from typing import List


class PengRobinsonEOS(EOSBase):

    @property
    def eos_u(self) -> float:
        return 2.0

    @property
    def eos_w(self) -> float:
        return -1.0

    @property
    def m_arr(self) -> np.ndarray:
        return 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
    
    def alpha(self, comp: Component) -> float:
        Tr = self.T / comp.Tc
        m = 0.37464 + 1.54226 * comp.omega - 0.26992 * comp.omega**2
        return (1 + m * (1 - math.sqrt(Tr)))**2

    def dalpha_dT(self, comp: Component) -> float:
        Tr = self.T / comp.Tc
        sqrt_Tr = math.sqrt(Tr)
        m = 0.37464 + 1.54226 * comp.omega - 0.26992 * comp.omega**2
        f = 1 + m * (1 - sqrt_Tr)
        return 2.0 * f * (-m / (2.0 * sqrt_Tr * comp.Tc))

    def alpha_array(self) -> np.ndarray:
        Tr = self.T / self.Tc
        sqrt_Tr = np.sqrt(Tr)
        return (1.0 + self.m_arr * (1.0 - sqrt_Tr))**2

    def dalpha_dT_array(self) -> np.ndarray:
        Tr = self.T / self.Tc
        sqrt_Tr = np.sqrt(Tr)
        f = 1.0 + self.m_arr * (1.0 - sqrt_Tr)
        return 2.0 * f * (-self.m_arr / (2.0 * sqrt_Tr * self.Tc))

    def a_formula(self, comp: Component) -> float:
        return 0.45724 * self.R**2 * comp.Tc**2 / comp.Pc

    def b_formula(self, comp: Component) -> float:
        return 0.07780 * self.R * comp.Tc / comp.Pc

    def get_coefficients(self, A: float, B: float) -> List[float]:
        return [1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]