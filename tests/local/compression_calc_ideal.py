from dataclasses import dataclass
from typing import Dict
import math

# ----------------------------------------------------------------------
# 1. Baza podataka za pojedine plinove (jako pojednostavljeno!)
#    M  [kg/kmol]
#    cp [kJ/(kg·K)] – približne vrijednosti
# ----------------------------------------------------------------------

SPECIES_DB = {
    "CO2": {"M": 44.01,   "cp": 0.846},
    "O2":  {"M": 31.999,  "cp": 0.918},
    "N2":  {"M": 28.0134, "cp": 1.04},
    "Ar":  {"M": 39.948,  "cp": 0.5203},  # kJ/(kg·K) @ 300 K :contentReference[oaicite:0]{index=0}
    "H2O": {"M": 18.015,  "cp": 1.87},    # “steam” ideal gas @ 300 K (tablica A-2) :contentReference[oaicite:1]{index=1}
    "NO":  {"M": 30.006,  "cp": 0.995},   # ~20 °C, 1 atm :contentReference[oaicite:2]{index=2}
    "SO2": {"M": 64.066,  "cp": 0.64},    # ~20 °C, 1 atm :contentReference[oaicite:3]{index=3}
    "SO3": {"M": 80.063,  "cp": 0.632},   # iz NIST Cp(298 K)=50.63 J/molK → 0.632 kJ/kgK :contentReference[oaicite:4]{index=4}
    "CO":  {"M": 28.01,   "cp": 1.04},    # @ 300 K, ista Wark/Kyle tablica kao za zrak itd. :contentReference[oaicite:5]{index=5}
    "H2S": {"M": 34.076,  "cp": 0.243}    # ~20 °C, 1 atm 
}

R_UNIV = 8.314  # kJ/(kmol·K)


# ----------------------------------------------------------------------
# 2. Klasa za smjesu idealnog plina
#    composition: dict("komponenta": molni_udio)
# ----------------------------------------------------------------------

@dataclass
class GasMixture:
    composition: Dict[str, float]  # molni udjeli

    def molar_mass(self) -> float:
        """Mješovita molarna masa [kg/kmol]."""
        M_mix = 0.0
        for comp, z in self.composition.items():
            sp = SPECIES_DB[comp]
            M_mix += z * sp["M"]
        return M_mix

    def cp(self) -> float:
        """
        Maseni cp smjese [kJ/(kg·K)].
        Kombiniramo cp po molnim udjelima (aproksimacija).
        """
        cp_molar = 0.0  # [kJ/(kmol·K)]
        for comp, z in self.composition.items():
            sp = SPECIES_DB[comp]
            cp_molar += z * sp["cp"] * sp["M"]  # cp_mass * M → cp_molar
        M_mix = self.molar_mass()
        cp_mass = cp_molar / M_mix            # [kJ/(kg·K)]
        return cp_mass

    def R(self) -> float:
        """Specifična plinska konstanta smjese [kJ/(kg·K)]."""
        M_mix = self.molar_mass()             # [kg/kmol]
        return R_UNIV / M_mix                 # [kJ/(kg·K)]

    def k(self) -> float:
        """Omjer toplinskih kapaciteta k = cp/cv [-]."""
        cp = self.cp()
        R = self.R()
        cv = cp - R
        return cp / cv


# ----------------------------------------------------------------------
# 3. Adijabatska (idealno izentropska + opcionalno stvarni η_s) kompresija
# ----------------------------------------------------------------------

def adiabatic_compression(
    p1: float,
    T1: float,
    p2: float,
    mixture: GasMixture,
    m_dot: float,
    eta_s: float = 1.0,
) -> dict:
    """
    Adijabatska kompresija idealnog plina s konstantnim cp.

    p1, p2  – tlakovi (u istim jedinicama, npr. bar ili Pa)
    T1      – ulazna temperatura [K]
    m_dot   – maseni protok [kg/s]
    eta_s   – izentropski stupanj djelovanja kompresora (1.0 = idealno)

    Vraća dict sa glavnim rezultatima.
    """
    cp = mixture.cp()  # [kJ/(kg·K)]
    k = mixture.k()

    # 1) Idealna izentropska adijabatska kompresija
    p_ratio = p2 / p1
    T2s = T1 * (p_ratio) ** ((k - 1.0) / k)      # [K]
    h1 = cp * T1                                 # [kJ/kg]
    h2s = cp * T2s                               # [kJ/kg]
    w_s = h2s - h1                               # [kJ/kg]

    # 2) Stvarni slučaj (ako eta_s < 1)
    if eta_s <= 0 or eta_s > 1.0:
        raise ValueError("eta_s treba biti u intervalu (0, 1].")

    w_actual = w_s / eta_s                      # [kJ/kg]
    h2 = h1 + w_actual                          # [kJ/kg]
    T2 = h2 / cp                                # [K] (jer h = cp T)

    # 3) Snaga
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


# ----------------------------------------------------------------------
# 4. Politropska kompresija idealnog plina
# ----------------------------------------------------------------------

def polytropic_compression(
    p1: float,
    T1: float,
    p2: float,
    mixture: GasMixture,
    m_dot: float,
    n: float,
) -> dict:
    """
    Politropska kompresija idealnog plina s eksponentom n (p v^n = const).

    p1, p2  – tlakovi (u istim jedinicama, npr. bar ili Pa)
    T1      – ulazna temperatura [K]
    m_dot   – maseni protok [kg/s]
    n       – politropski eksponent (između 1 i k)

    Vraća dict sa: T2, h1, h2, w, P_MW, p_ratio.
    """
    if n <= 1.0:
        raise ValueError("Politropski eksponent n mora biti > 1.")

    R = mixture.R()   # [kJ/(kg·K)]
    cp = mixture.cp() # [kJ/(kg·K)]

    p_ratio = p2 / p1

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


# ----------------------------------------------------------------------
# 5. Lijepi ispisi rezultata (na hrvatskom, s jedinicama)
# ----------------------------------------------------------------------

def print_adiabatic_results(p1, T1, p2, m_dot, res):
    print("==============================================")
    print(" ADIJABATSKA KOMPRESIJA IDEALNOG PLINA")
    print("==============================================")
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


# ----------------------------------------------------------------------
# 6. Primjer korištenja
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # smjesa npr. zraka ~ 79% N2, 21% O2
    #gass = GasMixture(composition={"N2": 0.79, "O2": 0.21})
    gass = GasMixture(composition={"CO2":0.85,"O2":0.0469, "N2":0.058, "Ar":0.0447, "H2O":0.0001, "NO":0.0001, "SO2":0.00005, "SO3":0.00005, "CO":0.00005, "H2S":0.00005})

    p1 = 1.0     # npr. bar
    T1 = 300.0   # K
    p2 = 100.0    # bar
    m_dot = 10.0 # kg/s

    # Adijabatski (izentropski i stvarni s eta_s=0.8)
    res_adiabatic = adiabatic_compression(p1, T1, p2, gass, m_dot, eta_s=0.8)
    print_adiabatic_results(p1, T1, p2, m_dot, res_adiabatic)

    # Politropski (n=1.3)
    n_poly = 1.2
    res_poly = polytropic_compression(p1, T1, p2, gass, m_dot, n=n_poly)
    print_polytropic_results(p1, T1, p2, m_dot, n_poly, res_poly)
