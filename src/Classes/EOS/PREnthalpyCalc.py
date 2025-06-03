import math
from typing import List
from src.Classes.Component import Component
from src.Classes.EOS.PengRobinsonEOS import PengRobinsonEOS
from src.Classes.EOS.RachfordRice import RachfordRice

class PREnthalpyCalc:

    R = 8.314

    def __init__(self, components: List[Component]):
        self._components = components

    def alpha(self, T, Tc, omega):
        Tr = T / Tc
        m = 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2
        return (1 + m * (1 - math.sqrt(Tr))) ** 2

    def da_dT_component(self, T, Tc, Pc, omega):
        Tr = T / Tc
        m = 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2
        d_alpha_dT = (-m / (Tc * math.sqrt(Tr))) * (1 + m * (1 - math.sqrt(Tr)))
        a = 0.45724 * self.R ** 2 * Tc ** 2 / Pc
        return a * d_alpha_dT

    def a_i(self, T, Tc, Pc, omega):
        alpha_val = self.alpha(T, Tc, omega)
        a = 0.45724 * self.R ** 2 * Tc ** 2 / Pc
        return a * alpha_val

    def b_i(self, Tc, Pc):
        return 0.07780 * self.R * Tc / Pc

    def da_dT(self, x: List[float], T: float):
        da = 0.0
        for xi, comp in zip(x, self._components):
            da += xi ** 2 * self.da_dT_component(T, comp.Tc, comp.Pc, comp.omega)
        return da

    def a_mix(self, x: List[float], T: float):
        a_mix = 0.0
        for i, xi in enumerate(x):
            comp_i = self._components[i]
            ai = self.a_i(T, comp_i.Tc, comp_i.Pc, comp_i.omega)
            for j, xj in enumerate(x):
                comp_j = self._components[j]
                aj = self.a_i(T, comp_j.Tc, comp_j.Pc, comp_j.omega)
                a_mix += xi * xj * math.sqrt(ai * aj)
        return a_mix

    def b_mix(self, x: List[float]):
        b_mix = 0.0
        for xi, comp in zip(x, self._components):
            bi = self.b_i(comp.Tc, comp.Pc)
            b_mix += xi * bi
        return b_mix

    def ideal_enthalpy(self, T: float, comp: Component, T_ref: float = 298.15) -> float:
        A, B, C, D = comp.CpA, comp.CpB, comp.CpC, comp.CpD
        return (A * (T - T_ref) +
                0.5 * B * (T ** 2 - T_ref ** 2) +
                (1.0 / 3.0) * C * (T ** 3 - T_ref ** 3) +
                0.25 * D * (T ** 4 - T_ref ** 4))

    def residual_enthalpy(self, x: List[float], T: float, P: float, Z: float) -> float:
        a_m = self.a_mix(x, T)
        b_m = self.b_mix(x)
        da_dT_m = self.da_dT(x, T)

        A = a_m * P / (self.R ** 2 * T ** 2)
        B = b_m * P / (self.R * T)
        sqrt2 = math.sqrt(2)
        log_term = math.log((Z + (1 + sqrt2) * B) / (Z + (1 - sqrt2) * B))

        h_res = self.R * T * (Z - 1) - (T * da_dT_m - a_m) / (b_m * sqrt2) * log_term
        return h_res

    def get_enthalpy(self, x: List[float], T: float, P: float, Z: float) -> float:
        h_ideal = sum(xi * self.ideal_enthalpy(T, comp) for xi, comp in zip(x, self._components))
        h_res = self.residual_enthalpy(x, T, P, Z)
        return h_ideal + h_res

    def get_total_enthalpy(self, x: List[float], y: List[float], Z_l: float, Z_v: float,
                            V: float, T: float, P: float) -> float:
        if V <= 0.0:
            return self.get_enthalpy(x, T, P, Z_l)
        elif V >= 1.0:
            return self.get_enthalpy(y, T, P, Z_v)
        else:
            H_liq = self.get_enthalpy(x, T, P, Z_l)
            H_vap = self.get_enthalpy(y, T, P, Z_v)
            return V * H_vap + (1 - V) * H_liq

    def compute_enthalpy_from_z(self, z: List[float], T: float, P: float) -> float:
        eos = PengRobinsonEOS(self._components, T, P)
        V, x, y, m, i = RachfordRice.solve(z, PengRobinsonEOS, self._components, T, P)
        Z_l, Z_v = eos.get_Z_factors(x, y)
        return self.get_total_enthalpy(x, y, Z_l, Z_v, V, T, P)
