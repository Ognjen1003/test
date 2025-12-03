import math
from .EOS import EOSBase
import src.Models.Component as Component
from typing import List
import numpy as np
import src.EnumsClasses.MethodsAndTypes as MT


class PengRobinsonEOS(EOSBase):
    
    def alpha(self, comp: Component) -> float:
        Tr = self.T / comp.Tc
        m = 0.37464 + 1.54226 * comp.omega - 0.26992 * comp.omega**2
        return (1 + m * (1 - np.sqrt(Tr)))**2

    def a_formula(self, comp: Component) -> float:
        return 0.45724 * self.R**2 * comp.Tc**2 / comp.Pc

    def b_formula(self, comp: Component) -> float:
        return 0.07780 * self.R * comp.Tc / comp.Pc

    def get_coefficients(self, A: float, B: float) -> List[float]:
        return [1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]

    def h_residual(self, p: float, T: float) -> float:
        """
        Peng–Robinson enthalpy departure h^R [kJ/kg]
        Formula koristi Z, A, B i d(alpha)/dT.
        """
        
        T_old, P_old = self.T, self.P
        self.T, self.P = T, p

        # x iz components
        # mix parametri
        x = np.array([c.fraction for c in self.components], dtype=float)
        a_mix, b_mix, A, B = self.get_mixture_parameters(x)

        # Z (vapor root!)
        Z = self.calc_Z_factor(A, B, MT.Phase.VAPOR)

        # Peng–Robinson parameters
        sqrt2 = math.sqrt(2)
        m = [0.37464 + 1.54226*c.omega - 0.26992*c.omega**2 for c in self.components]

        # alpha_i and its T derivative
        dalpha_dT = []
        alpha = []
        for i, c in enumerate(self.components):
            Tr = T / c.Tc
            m_i = m[i]
            sqrt_Tr = math.sqrt(Tr)
            a_i_alpha = self.a_formula(c) * (1 + m_i*(1 - sqrt_Tr))**2
            alpha.append(a_i_alpha)
            # derivative wrt T:
            dalpha_i = self.a_formula(c) * 2*(1 + m_i*(1 - sqrt_Tr)) * (-m_i/(2*sqrt_Tr)) * (1/c.Tc)
            dalpha_dT.append(dalpha_i)

        # mixture derivative d(a_mix)/dT
        sqrt_ai = np.sqrt(np.array(alpha))
        u = x * sqrt_ai

        if self.k_ij is None:
            da_mix_dT = 2 * np.sum(u) * np.sum(0.5 * x / sqrt_ai * np.array(dalpha_dT))
        else:
            K = np.asarray(self.k_ij)
            t = (1 - K) @ u
            du_dT = x * (0.5/ sqrt_ai) * np.array(dalpha_dT)
            da_mix_dT = 2 * np.dot(du_dT, t)

        # actual PR formula:
        term1 = Z - 1
        term2 = -(A / (2*B*sqrt2)) * (da_mix_dT / a_mix) * math.log((Z + (1+sqrt2)*B)/(Z + (1-sqrt2)*B))

        hR = self.R * T * (term1 + term2)  # [per mole]

        # convert to per kg:
        M_mix = sum(ci.fraction * ci.Mw for ci in self.components)
        hR_mass = hR / M_mix

        self.T, self.P = T_old, P_old

        return hR_mass
    
    def s_residual(self, p: float, T: float) -> float:
        """
        Peng–Robinson entropy departure s^R [kJ/(kg·K)]
        """
        T_old, P_old = self.T, self.P
        self.T, self.P = T, p

        x = np.array([c.fraction for c in self.components], dtype=float)
        
        a_mix, b_mix, A, B = self.get_mixture_parameters(x)
        Z = self.calc_Z_factor(A, B, MT.Phase.VAPOR)
        sqrt2 = math.sqrt(2)

        term1 = math.log(Z - B)
        term2 = -(A / (2*B*sqrt2)) * math.log((Z + (1+sqrt2)*B)/(Z + (1-sqrt2)*B))

        sR = self.R * (term1 + term2)   # per mole

        M_mix = sum(ci.fraction * ci.Mw for ci in self.components)
        sR_mass = sR / M_mix

        self.T, self.P = T_old, P_old

        return sR_mass


