from src.Classes.EOS.PengRobinsonEOS import PengRobinsonEOS
from src.Classes.EOS.RachfordRice import RachfordRice
from src.Models.Component import Component
import src.EnumsClasses.MethodsAndTypes as MT
from typing import List, Dict
from dataclasses import dataclass



@dataclass
class State:
    p: float
    T: float
    h: float
    s: float
    v: float
    z: Dict[str, float]

@dataclass
class GasMixture:
    composition: List[Component]  # molni udjeli

    def molar_mass(self) -> float:
        """Mješovita molarna masa [kg/kmol]."""
        M_mix = 0.0
        for comp in self.composition:
            M_mix += comp.fraction * comp.Mw
        return M_mix

    def cp(self) -> float:
        """
        Maseni cp smjese [kJ/(kg·K)].
        Kombiniramo cp po molnim udjelima (aproksimacija).
        """
        cp_molar = 0.0  # [kJ/(kmol·K)]
        for comp in self.composition:
            cp_molar += comp.fraction * comp.Cp * comp.Mw  # cp_mass * M → cp_molar
        M_mix = self.molar_mass()
        cp_mass = cp_molar / M_mix            # [kJ/(kg·K)]
        return cp_mass

    def R(self) -> float:
        """Specifična plinska konstanta smjese [kJ/(kg·K)]."""
        M_mix = self.molar_mass()             # [kg/kmol]
        return MT.CONSTANTS.R / M_mix                 # [kJ/(kg·K)]

    def k(self) -> float:
        """Omjer toplinskih kapaciteta k = cp/cv [-]."""
        cp = self.cp()
        R = self.R()
        cv = cp - R
        return cp / cv

@dataclass
class CompressorThermo:

    class Ideal:

        @staticmethod
        def calc_ideal_gas_sanity_check(P1: int, P2: int, T1: int, mass_flow: float, isentropic_efficiency: float, polytropic_exponent: float, gass_data: List[Component]):
            
            gass = GasMixture(gass_data)
            cp = gass.cp()  # [kJ/(kg·K)]
            k = gass.k()
            R = gass.R()
            p_ratio = P2 / P1

            # Adijabatski (izentropski i stvarni s eta_s)
            res_adiabatic = CompressorThermo.Ideal.adiabatic_compression(P1, T1, P2, cp, k, p_ratio, mass_flow, isentropic_efficiency)
            #print_adiabatic_results(P1, T1, P2, mass_flow, res_adiabatic, isentropic_efficiency)

            # Politropski (n od 1-k)
            res_poly = CompressorThermo.Ideal.polytropic_compression(P1, T1, P2, cp, R, p_ratio, mass_flow, polytropic_exponent)
            #print_polytropic_results(P1, T1, P2, mass_flow, n_poly, res_poly)
            
            
            return {
                "adiabatic_ideal": res_adiabatic,
                "polytropic_ideal": res_poly
            }


        # ----------------------------------------------------------------------
        # 3. Adijabatska (idealno izentropska + opcionalno stvarni η_s) kompresija
        # ----------------------------------------------------------------------
        @staticmethod
        def adiabatic_compression( p1: float, T1: float, p2: float, cp: float, k: float, p_ratio: float, m_dot: float, eta_s: float = 1.0 ) -> dict:
            """
            Adijabatska kompresija idealnog plina s konstantnim cp.
            p1, p2  tlakovi (u istim jedinicama, npr. bar ili Pa)
            T1      ulazna temperatura [K]
            m_dot   maseni protok [kg/s]
            eta_s   izentropski stupanj djelovanja kompresora (1.0 = idealno)
            """
  
            T2s = T1 * (p_ratio) ** ((k - 1.0) / k)      # [K]
            h1 = cp * T1                                 # [kJ/kg]
            h2s = cp * T2s                               # [kJ/kg]
            w_s = h2s - h1                               # [kJ/kg]

            if eta_s <= 0 or eta_s > 1.0:
                raise ValueError("eta_s treba biti u intervalu (0, 1].")

            w_actual = w_s / eta_s                      # [kJ/kg]
            h2 = h1 + w_actual                          # [kJ/kg]
            T2 = h2 / cp                                # [K] (jer h = cp T)

            # Snaga
            P_s_kW = m_dot * w_s                        # [kW]
            P_kW = m_dot * w_actual                     # [kW]
            P_s_MW = P_s_kW / 1e3                       # [MW]
            P_MW = P_kW / 1e3                           # [MW]

            return {
                "T2s": T2s,
                "T2": T2,
                "h1": h1,
                "h2s": h2s,
                "h2": h2,
                "w_s": w_s,
                "w_actual": w_actual,
                "P_s_MW": P_s_MW,
                "P_MW": P_MW,
                "p_ratio": p_ratio,
            }

        @staticmethod
        def polytropic_compression( p1: float, T1: float, p2: float, cp: float, R: float, p_ratio: float, m_dot: float, n: float) -> dict:
            """
            Politropska kompresija idealnog plina s eksponentom n (p v^n = const).

            p1, p2  tlakovi (u istim jedinicama, npr. bar ili Pa)
            T1      ulazna temperatura [K]
            m_dot   maseni protok [kg/s]
            n       politropski eksponent (između 1 i k)

            Vraća dict sa: T2, h1, h2, w, P_MW, p_ratio.
            """
            if n <= 1.0:
                raise ValueError("Politropski eksponent n mora biti > 1.")

            # 1) T2 iz politropske relacije za idealni plin
            T2 = T1 * (p_ratio) ** ((n - 1.0) / n)      # [K]

            # 2) Rad politropske kompresije za idealni plin
            #    w = (n/(n-1)) * R * T1 * [ (p2/p1)^((n-1)/n) - 1 ]
            w = (n / (n - 1.0)) * R * T1 * ((p_ratio) ** ((n - 1.0) / n) - 1.0)  # [kJ/kg]

            # 3) Enthalpije po idealnom plinu (h = cp T)
            h1 = cp * T1                                # [kJ/kg]
            h2 = cp * T2                                # [kJ/kg]

            # 4) Snaga
            P_kW = m_dot * w                            # [kW]
            P_MW = P_kW / 1e3                           # [MW]

            return {
                "T2": T2,
                "h1": h1,
                "h2": h2,
                "w": w,
                "P_MW": P_MW,
                "p_ratio": p_ratio,
            }

        @staticmethod
        def print_adiabatic_results(p1, T1, p2, m_dot, eta_s, res ):
            print("==============================================")
            print(" ADIJABATSKA KOMPRESIJA IDEALNOG PLINA")
            print("==============================================")
            print(f"eta_s {eta_s}")
            print(f"Ulazni tlak p1:        {p1:.3f} [isto kao p2]")
            print(f"Izlazni tlak p2:       {p2:.3f} [isto kao p1]")
            print(f"Omjer tlakova p2/p1:   {res['p_ratio']:.3f} [-]\n")

            print(f"Ulazna temperatura T1: {T1:.2f} K ({T1 - 273.15:.2f} °C)")
            print(f"Izlaz T2s (idealno, izentropski): {res['T2s']:.2f} K ({res['T2s'] - 273.15:.2f} °C)")
            print(f"Izlaz T2 (stvarni proces):        {res['T2']:.2f} K ({res['T2'] - 273.15:.2f} °C)\n")

            print("Specifične entalpije (referentna nula proizvoljna):")
            print(f"  h1  (ulaz):                      {res['h1']:.3f} kJ/kg")
            print(f"  h2s (idealni izentropski izlaz): {res['h2s']:.3f} kJ/kg")
            print(f"  h2  (stvarni izlaz):             {res['h2']:.3f} kJ/kg\n")

            print("Specifični rad kompresije:")
            print(f"  w_s      (idealni, izentropski): {res['w_s']:.3f} kJ/kg")
            print(f"  w_actual (stvarni proces):       {res['w_actual']:.3f} kJ/kg\n")

            print(f"Maseni protok:        {m_dot:.3f} kg/s")
            print(f"P_s (idealna snaga):  {res['P_s_MW']:.4f} MW")
            print(f"P   (stvarna snaga):  {res['P_MW']:.4f} MW")
            print("==============================================\n")

        @staticmethod
        def print_polytropic_results(p1, T1, p2, m_dot, n, res):
            print("==============================================")
            print(" POLITROPSKA KOMPRESIJA IDEALNOG PLINA")
            print("==============================================")
            print(f"Ulazni tlak p1:        {p1:.3f} [isto kao p2]")
            print(f"Izlazni tlak p2:       {p2:.3f} [isto kao p1]")
            print(f"Omjer tlakova p2/p1:   {res['p_ratio']:.3f} [-]\n")

            print(f"Politropski eksponent n: {n:.3f} [-]\n")

            print(f"Ulazna temperatura T1: {T1:.2f} K ({T1 - 273.15:.2f} °C)")
            print(f"Izlazna temperatura T2: {res['T2']:.2f} K ({res['T2'] - 273.15:.2f} °C)\n")

            print("Specifične entalpije (idealni plin, h = cp·T):")
            print(f"  h1 (ulaz):   {res['h1']:.3f} kJ/kg")
            print(f"  h2 (izlaz):  {res['h2']:.3f} kJ/kg\n")

            print("Specifični rad politropske kompresije:")
            print(f"  w:           {res['w']:.3f} kJ/kg\n")

            print(f"Maseni protok: {m_dot:.3f} kg/s")
            print(f"Snaga P:       {res['P_MW']:.4f} MW")
            print("==============================================\n")