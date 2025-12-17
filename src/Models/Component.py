import src.EnumsClasses.MethodsAndTypes as MT
from pydantic import BaseModel
import math

class Component(BaseModel):
    name: str
    formula: str
    Mw: float       # [kg/kmol]
    Tc: float
    Pc: float
    omega: float
    fraction: float
    Cp: float | None = None 

    # NASA 7-coefficient polynomials
    nasa_low: list[float] | None = None    # [a1..a7]
    nasa_high: list[float] | None = None   # [a1..a7]
    nasa_mid: float | None = None          # granica low→high (najčešće 1000 K)

    def _get_nasa9_coeffs(self, T: float):
        """
        Vrati (a1..a7, b1, b2) za zadanu temperaturu T
        iz low ili high seta, ovisno o T_mid.
        """

        if self.nasa_mid is None:
            a1, a2, a3, a4, a5, a6, a7, b1, b2 = self.nasa_low

        if T <= self.nasa_mid:
            a1, a2, a3, a4, a5, a6, a7, b1, b2 = self.nasa_low
        else:
            a1, a2, a3, a4, a5, a6, a7, b1, b2 = self.nasa_high
        return a1, a2, a3, a4, a5, a6, a7, b1, b2

    def cp_ideal_molar(self, T: float) -> float:
        """
        cp°(T) po MOLU [kJ/(kmol·K)] prema NASA-9.
        """
        a1, a2, a3, a4, a5, a6, a7, b1, b2 = self._get_nasa9_coeffs(T)
        R_univ = 8.314462618  # kJ/(kmol·K)

        cpT_over_R = (
            a1 * T**-2
            + a2 * T**-1
            + a3
            + a4 * T
            + a5 * T**2
            + a6 * T**3
            + a7 * T**4
        )
        cp_over_R = cpT_over_R / T
        return cp_over_R * R_univ

    def cp_ideal_mass(self, T: float) -> float:
        """
        cp°(T) po MASI [kJ/(kg·K)].
        """
        return self.cp_ideal_molar(T) / self.Mw

    def h_ideal_molar(self, T: float) -> float:
        """
        h°(T) po MOLU [kJ/kmol] prema NASA-9.
        """
        a1, a2, a3, a4, a5, a6, a7, b1, b2 = self._get_nasa9_coeffs(T)
        R_univ = 8.314462618  # kJ/(kmol·K)

        h_over_RT = (
            -a1 * T**-2
            + a2 * math.log(T) / T
            + a3
            + a4 * T / 2.0
            + a5 * T**2 / 3.0
            + a6 * T**3 / 4.0
            + a7 * T**4 / 5.0
            + b1 / T
        )
        return h_over_RT * R_univ * T  # kJ/kmol

    def h_ideal_mass(self, T: float) -> float:
        """
        h°(T) po MASI [kJ/kg].
        """
        return self.h_ideal_molar(T) / self.Mw

    def s_ideal_molar(self, T: float, p: float, pref: float = 1e5) -> float:
        """
        s°(T,p) po MOLU [kJ/(kmol·K)] prema NASA-9.
        """
        a1, a2, a3, a4, a5, a6, a7, b1, b2 = self._get_nasa9_coeffs(T)

        s_over_R_T = (
            -a1 * T**-2 / 2.0
            - a2 * T**-1
            + a3 * math.log(T)
            + a4 * T
            + a5 * T**2 / 2.0
            + a6 * T**3 / 3.0
            + a7 * T**4 / 4.0
            + b2
        )

        s_over_R = s_over_R_T - math.log(p / pref)
        return s_over_R * MT.CONSTANTS.R  # kJ/(kmol·K)

    def s_ideal_mass(self, T: float, p: float, pref: float = 1e5) -> float:
        """
        s°(T,p) po MASI [kJ/(kg·K)].
        """
        return self.s_ideal_molar(T, p, pref) / self.Mw

