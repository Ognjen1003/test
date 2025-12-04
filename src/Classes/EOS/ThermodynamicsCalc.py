from src.Classes.EOS.PengRobinsonEOS import PengRobinsonEOS
from src.Classes.EOS.RachfordRice import RachfordRice
from src.Models.Component import Component
import src.EnumsClasses.MethodsAndTypes as MT
from typing import List, Dict, Callable
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
                "n": n,
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

    class Real:

        @staticmethod
        def calc_real_gas_thermo(P1: int, P2: int, T1: int, m_dot: float, isentropic_efficiency: float, polytropic_efficiency: float, gass_data: List[Component]):
            
            eos = PengRobinsonEOS(gass_data, T1, P1)

            res_ad = CompressorThermo.Real.adiabatic_compression_real(P1, T1, P2, m_dot, eos, isentropic_efficiency)
            res_poly = CompressorThermo.Real.polytropic_compression_real(P1, T1, P2, m_dot, eos, polytropic_efficiency)

            #CompressorThermo.Real.print_adiabatic_results(P1, T1, P2, m_dot, isentropic_efficiency, res_ad )
            #CompressorThermo.Real.print_polytropic_results(P1, T1, P2, m_dot, polytropic_efficiency, res_poly )
            
            return {
                "adiabatic_ideal": res_ad,
                "polytropic_ideal": res_poly
            }


        @staticmethod
        def build_state(p: float, T: float, eos: PengRobinsonEOS) -> State:
            """Iz tlaka, temperature i kompozicije izgradi kompletno stanje."""
            h = eos.h(p, T)
            s = eos.s(p, T)
            return State(p=p, T=T, h=h, s=s, v=None, z=eos.components)
        
        # ================================================================
        #     Jednostavan 1D root-finder po temperaturi
        #    (brutalno jednostavno, ali dovoljno za strukturu)
        # ================================================================
        @staticmethod
        def find_T_for_target(
            func_T: Callable[[float], float],
            T_min: float,
            T_max: float,
            tol: float = 1e-4,
            max_iter: int = 50,
        ) -> float:
            """
            Bisection za rješavanje func_T(T) = 0 u [T_min, T_max].
            Koristi se za:
            - s(p, T) - s_target = 0
            - h(p, T) - h_target = 0
            Pretpostavka: funkcija mijenja predznak na krajevima intervala.
            """
            f_min = func_T(T_min)
            f_max = func_T(T_max)

            if f_min * f_max > 0:
                raise ValueError("Root-finder: funkcija nema promjenu predznaka u zadanom intervalu.")

            for _ in range(max_iter):
                T_mid = 0.5 * (T_min + T_max)
                f_mid = func_T(T_mid)

                if abs(f_mid) < tol:
                    return T_mid

                if f_min * f_mid < 0:
                    T_max = T_mid
                    f_max = f_mid
                else:
                    T_min = T_mid
                    f_min = f_mid

            # Ako nismo konvergirali, svejedno vraćamo sredinu
            return 0.5 * (T_min + T_max)


        # Helperi specifični za EOS:
        @staticmethod
        def solve_T_at_const_s(p: float, s_target: float, eos: PengRobinsonEOS,
                            T_min: float, T_max: float) -> float:
            """Nađi T tako da s(p, T, z) = s_target."""
            def f(T):
                return eos.s(p, T) - s_target
            return CompressorThermo.Real.find_T_for_target(f, T_min, T_max)

        @staticmethod
        def solve_T_at_const_h(p: float, h_target: float, eos: PengRobinsonEOS,
                            T_min: float, T_max: float) -> float:
            """Nađi T tako da h(p, T, z) = h_target."""
            def f(T):
                return eos.h(p, T) - h_target
            return CompressorThermo.Real.find_T_for_target(f, T_min, T_max)
        
        @staticmethod
        def adiabatic_compression_real(
            p1: float,
            T1: float,
            p2: float,
            m_dot: float,
            eos: PengRobinsonEOS,
            eta_s: float = 1.0,
            T_bracket: tuple = (250.0, 1500.0),
        ) -> dict:
            """
            Adijabatska kompresija realne smjese (PR/GERG EOS).

            p1, p2   – tlak (npr. bar ili Pa, ali konzistentno u cijelom kodu)
            T1       – ulazna temperatura [K]
            z        – kompozicija smjese (molni ili maseni udjeli, konzistentno s EOS-om)
            m_dot    – maseni protok [kg/s]
            eos      – objekt koji implementira EoSModel (tvoj PR kod)
            eta_s    – izentropski stupanj djelovanja kompresora (1.0 = idealno)
            T_bracket – interval [T_min, T_max] unutar kojeg tražimo rješenja za T [K].

            Vraća dict s ključnim veličinama.
            """

            if not (0 < eta_s <= 1.0):
                raise ValueError("eta_s treba biti u intervalu (0, 1].")

            # 1) Početno stanje
            st1 = CompressorThermo.Real.build_state(p1, T1, eos)

            # 2) Idealno izentropsko izlazno stanje: s2s = s1, p2 zadano
            s1 = st1.s
            T2s = CompressorThermo.Real.solve_T_at_const_s(
                p=p2,
                s_target=s1,
                eos=eos,
                T_min=T_bracket[0],
                T_max=T_bracket[1],
            )
            st2s = CompressorThermo.Real.build_state(p2, T2s, eos)

            # 3) Idealni izentropski specifični rad
            w_s = st2s.h - st1.h   # [kJ/kg]

            # 4) Stvarni izlaz s eta_s
            #    w_actual = w_s / eta_s → h2 = h1 + w_actual
            w_actual = w_s / eta_s
            h2 = st1.h + w_actual

            T2 = CompressorThermo.Real.solve_T_at_const_h(
                p=p2,
                h_target=h2,
                eos=eos,
                T_min=T_bracket[0],
                T_max=T_bracket[1],
            )
            st2 = CompressorThermo.Real.build_state(p2, T2, eos)

            # 5) Snaga
            P_s_kW = m_dot * w_s          # [kW]
            P_kW = m_dot * w_actual       # [kW]
            P_s_MW = P_s_kW / 1e3
            P_MW = P_kW / 1e3

            return {
                "state1": st1,
                "state2s": st2s,
                "state2": st2,
                "w_s": w_s,
                "w_actual": w_actual,
                "P_s_MW": P_s_MW,
                "P_MW": P_MW,
                "p_ratio": p2 / p1,
                "eta_s": eta_s,
            }

        @staticmethod
        def polytropic_compression_real(
            p1: float,
            T1: float,
            p2: float,
            m_dot: float,
            eos: PengRobinsonEOS,
            eta_p: float,
            n_steps: int = 50,
            T_bracket: tuple = (250.0, 1500.0),
        ) -> dict:
            """
            Politropska kompresija realne smjese korištenjem diskretizacije po tlaku.

            Ideja:
            - Podijelimo raspon tlaka [p1, p2] na n_steps
            - U svakom koraku:
                * idealni izentropski skok pri Δp
                * korekcija na stvarni politropski skok preko η_p
            - w_total = h2 - h1 (zbroj Δh_actual)

            Parametri:
            p1, p2   – početni i krajnji tlak
            T1       – početna temperatura [K]
            z        – kompozicija smjese
            m_dot    – maseni protok [kg/s]
            eos      – PR/GERG model
            eta_p    – politropska učinkovitost kompresora (0–1)
            n_steps  – broj diskretnih koraka po tlaku
            T_bracket – interval za traženje T po koracima

            Vraća dict s ukupnim radom, snagom i početnim/završnim stanjem.
            """

            if not (0 < eta_p <= 1.0):
                raise ValueError("eta_p treba biti u intervalu (0, 1].")

            # 1) Početno stanje
            st = CompressorThermo.Real.build_state(p1, T1, eos)
            st1 = st
            p_start = p1
            p_end = p2

            # Pretpostavimo p2 > p1 (kompresija). Ako nije, može se dodati provjera.
            dp = (p_end - p_start) / n_steps

            for i in range(n_steps):
                p_old = p_start + i * dp
                p_new = p_old + dp

                # 2) Idealni izentropski skok u ovom malom koraku:
                s_old = st.s
                T_iso = CompressorThermo.Real.solve_T_at_const_s(
                    p=p_new,
                    s_target=s_old,
                    eos=eos,
                    T_min=T_bracket[0],
                    T_max=T_bracket[1],
                )
                st_iso = CompressorThermo.Real.build_state(p_new, T_iso, eos)

                Δh_s = st_iso.h - st.h    # idealni Δh u tom koraku

                # 3) Politropska korekcija sa η_p:
                #    Δh_actual = Δh_s / η_p
                Δh_actual = Δh_s / eta_p
                h_new = st.h + Δh_actual

                # 4) Nađi T_new iz uvjeta h(p_new, T_new, z) = h_new
                T_new = CompressorThermo.Real.solve_T_at_const_h(
                    p=p_new,
                    h_target=h_new,
                    eos=eos,
                    T_min=T_bracket[0],
                    T_max=T_bracket[1],
                )
                st = CompressorThermo.Real.build_state(p_new, T_new, eos)

            st2 = st

            # 5) Ukupni specifični rad i snaga
            w_total = st2.h - st1.h   # [kJ/kg]
            P_kW = m_dot * w_total
            P_MW = P_kW / 1e3

            return {
                "state1": st1,
                "state2": st2,
                "w": w_total,
                "P_MW": P_MW,
                "p_ratio": p2 / p1,
                "eta_p": eta_p,
                "n_steps": n_steps,
            }
        
        def print_adiabatic_results(p1, T1, p2, m_dot, eta_s, res):
            st1 = res["state1"]
            st2s = res["state2s"]
            st2 = res["state2"]

            print("\n====================================================")
            print(" ADIJABATSKA KOMPRESIJA REALNOG PLINA (PR EOS)")
            print("====================================================\n")

            print("► Ulazni parametri:")
            print(f"  - Ulazni tlak p1:                 {p1:.3f} bar")
            print(f"  - Izlazni tlak p2:                {p2:.3f} bar")
            print(f"  - Omjer tlakova p2/p1:            {res['p_ratio']:.3f} [-]")
            print(f"  - Ulazna temperatura T1:          {T1:.2f} K   ({T1 - 273.15:.2f} °C)")
            print(f"  - Maseni protok ṁ:               {m_dot:.3f} kg/s")
            print(f"  - Izentropski stupanj η_s:        {eta_s:.3f}\n")

            print("► Stanja izračunata EOS-om:")
            print(f"  - Izentropska izlazna T2s:        {st2s.T:.2f} K   ({st2s.T - 273.15:.2f} °C)")
            print(f"  - Stvarna izlazna temperatura T2: {st2.T:.2f} K   ({st2.T - 273.15:.2f} °C)\n")

            print("► Entalpije (apsolutna razina proizvoljna):")
            print(f"  - h1  (ulaz):                     {st1.h:.3f} kJ/kg")
            print(f"  - h2s (idealno izentropski):      {st2s.h:.3f} kJ/kg")
            print(f"  - h2  (stvarni proces):           {st2.h:.3f} kJ/kg\n")

            print("► Specifični rad kompresije:")
            print(f"  - w_s      (idealni, izentropski): {res['w_s']:.3f} kJ/kg")
            print(f"  - w_actual (stvarni):              {res['w_actual']:.3f} kJ/kg\n")

            print("► Snaga kompresora:")
            print(f"  - Idealna snaga P_s:               {res['P_s_MW']:.4f} MW")
            print(f"  - Stvarna snaga  P:                {res['P_MW']:.4f} MW")
            print("====================================================\n")



        def print_polytropic_results(p1, T1, p2, m_dot, eta_p, res):
            st1 = res["state1"]
            st2 = res["state2"]

            print("\n====================================================")
            print(" POLITROPSKA KOMPRESIJA REALNOG PLINA (PR EOS)")
            print("====================================================\n")

            print("► Ulazni parametri:")
            print(f"  - Ulazni tlak p1:                  {p1:.3f} bar")
            print(f"  - Izlazni tlak p2:                 {p2:.3f} bar")
            print(f"  - Omjer tlakova p2/p1:             {res['p_ratio']:.3f} [-]")
            print(f"  - Ulazna temperatura T1:           {T1:.2f} K ({T1 - 273.15:.2f} °C)")
            print(f"  - Maseni protok ṁ:                {m_dot:.3f} kg/s")
            print(f"  - Politropska učinkovitost η_p:    {eta_p:.3f}")
            print(f"  - Broj diskretizacijskih koraka:   {res['n_steps']}\n")

            print("► Rezultati stanja:")
            print(f"  - Izlazna temperatura T2:          {st2.T:.2f} K ({st2.T - 273.15:.2f} °C)")

            print("\n► Entalpije (prema EOS):")
            print(f"  - h1 (ulaz):                       {st1.h:.3f} kJ/kg")
            print(f"  - h2 (izlaz):                      {st2.h:.3f} kJ/kg")

            print("\n► Specifični rad politropske kompresije:")
            print(f"  - w_total:                         {res['w']:.3f} kJ/kg")

            print("\n► Snaga kompresora:")
            print(f"  - P:                               {res['P_MW']:.4f} MW")

            print("====================================================\n")