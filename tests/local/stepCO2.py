from typing import List, Dict
import src.Classes.EOS.RachfordRice as RR
from src.Classes.EOS import PengRobinsonEOS
from Models.Component import Component
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import pandas as pd

# ===== GENERALNE STVARI ZA SAD =====
R = 8.314462618  # J/(mol·K)
V_LOOP = 0.05  # m3
P_target_bar=40
T_target_C=40
P_init_bar=2.2
T_init_C=20
dosing_volume_liters=10

# === DEFINICIJA KOMPONENTI ===
components = [
    Component(name="CarbonDioxide", formula="CO2", Mw=44.01, Tc=304.7, Pc=73.866, omega=0.225, fraction=0.91, CpA=22.26, CpB=5.981e-3, CpC=-3.501e-6, CpD=7469.0),
    Component(name="Nitrogen", formula="N2", Mw=28.013, Tc=126.2, Pc=33.944, omega=0.04, fraction=0.04, CpA=28.90, CpB=-0.0001571, CpC=8.081e-7, CpD=0.0),
    Component(name="Argon", formula="Ar", Mw=39.948, Tc=150.8, Pc=48.63, omega=-0.004, fraction=0.003, CpA=20.85, CpB=0.0, CpC=0.0, CpD=0.0),
    Component(name="Water", formula="H2O", Mw=18.015, Tc=647.1, Pc=220.64, omega=0.344, fraction=0.002, CpA=32.24, CpB=0.0, CpC=0.0, CpD=0.0),
    Component(name="SulfurDioxide", formula="SO2", Mw=64.066, Tc=430.8, Pc=77.8, omega=0.251, fraction=0.045, CpA=40.0, CpB=0.0, CpC=0.0, CpD=0.0)
]

# === FUNKCIJA ZA NORMALIZACIJU MOLSKIH FRAKCIJA NA V_LOOP L UKUPNO (CoolProp) ===
def normalize_mols_to_volume(df_input, target_volume_liters):
    df = df_input.copy()
    df["V_target"] = pd.to_numeric(df["V_target"], errors='coerce')
    total_V = df["V_target"].sum()
    correction_factor = target_volume_liters / total_V
    df["Mol [mol]"] *= correction_factor
    df["V_target"] *= correction_factor
    df["V_CoolProp [L]"] *= correction_factor
    df["V_ideal [L]"] *= correction_factor
    df["Doses"] *= correction_factor
    df["Doses adjusted"] *= correction_factor
    return df

# === FUNKCIJA ZA PRORAČUN IDEALNOG I COOLPROP pretpostavke ===
def calculate_charging():
    P_target = P_target_bar * 1e5
    T_target = T_target_C + 273.15
    P_init = P_init_bar * 1e5
    T_init = T_init_C + 273.15

    print(f"================================= m3, Pa i K =================================")
    print(f"V_LOOP: {V_LOOP}")
    print(f"P_init: {P_init}, T_init: {T_init}")
    print(f"P_target: {P_target}, T_target: {T_target}")
    print(f"==============================================================================")

    n_total = P_target * V_LOOP / (R * T_target)
    print(f"n_total molovi u idealnom slucaju pv=nRT: {n_total}" )

    results = []
    for comp in components:

        n_i = comp.fraction * n_total

        if(comp.name == 'Water'):
            V_water = n_i * comp.Mw/1e6
            results.append({
            "Component": comp.name,
            "Mol [mol]": round(n_i, 4),
            "V_ideal [L]": round(V_water, 8),
            "V_CoolProp [L]": round(V_water, 8),
            "Doses": 0.0,
            "V_target": round(V_water, 8),
            "Doses adjusted": 0.0,
            })
            continue

        V_ideal = n_i * R * T_init / P_init

        try:
            rho_molar = PropsSI("Dmolar", "T", T_init, "P", P_init, comp.formula)  # mol/m3
            v_coolprop = 1 / rho_molar  # m3/mol
            V_CoolProp = n_i * v_coolprop         

            rho_molar_target = PropsSI("Dmolar", "T", T_target, "P", P_target, comp.formula)  # mol/m3
            v_coolprop_target = 1 / rho_molar_target  # m3/mol
            V_CoolProp_target = n_i * v_coolprop_target                       
        except:
            print(f"comp.name: {comp.name} nije u listi od CoolPropa, downgrade na Vm idealni")
            V_CoolProp = 22.4 * n_i
            V_CoolProp_target = 22.4 * n_i

        doses = V_ideal * 1000 / dosing_volume_liters
        doses_adjusted = V_CoolProp_target * 1000 / dosing_volume_liters
        results.append({
            "Component": comp.name,
            "Mol [mol]": round(n_i, 4),
            "V_ideal [L]": round(V_ideal * 1000, 4),
            "V_CoolProp [L]": round(V_CoolProp * 1000, 4),
            "Doses": round(doses, 4),
            "V_target": round(V_CoolProp_target * 1000, 4),
            "Doses adjusted": round(doses_adjusted, 4),
        })

    return pd.DataFrame(results)

def calculate_pressure_from_eos(df_input, T_K, V_m3):

    component_objs = []
    fractions = []
    for _, row in df_input.iterrows():
        comp = next((c for c in components if c.name == row["Component"]), None)
        if comp:
            c_new = Component(
                name=comp.name,
                formula=comp.formula,
                Mw=comp.Mw,
                Tc=comp.Tc,
                Pc=comp.Pc * 1e5,  #Pa
                omega=comp.omega,
                fraction=row["Mol [mol]"] / df_input["Mol [mol]"].sum(),
                CpA=comp.CpA,
                CpB=comp.CpB,
                CpC=comp.CpC,
                CpD=comp.CpD
            )
            component_objs.append(c_new)

    eos = PengRobinsonEOS(component_objs, T=T_K, P=1e5)  # dummy P, provjerit
    n_total = df_input["Mol [mol]"].sum()
    v_molar = V_m3 / n_total  # m3/mol
    # Phase = RachfordRice.return_phase()
    #if....
    #else...
    p_real = eos.get_pressure_with_Z(v_molar, phase='vapor')
    return p_real


if __name__ == "__main__":
    

   P_target = P_target_bar * 1e5
   T_target = T_target_C + 273.15
   P_init = P_init_bar * 1e5
   T_init = T_init_C + 273.15

   n_total = P_target * V_LOOP / (R * T_target)
   rho_molar = PropsSI("Dmolar", "T", T_init, "P", P_init, "CarbonDioxide")  # mol/m3
   v_coolprop = 1 / rho_molar  # m3/mol
   v_coolprop_litre_po_molu = v_coolprop * 1000 # l/mol

   n_CO2 = 0.91 * n_total  

   n_molova_iz_coolpropa = V_LOOP / v_coolprop

   print(f"n_total: {n_total}")    
   print(f"n_CO2: {n_CO2}")  
   print(f"v_coolprop_litre_po_molu: {v_coolprop_litre_po_molu}")  # molarni volumen pri tim uvjetima u litrama
   print(f"n_molova_iz_coolpropa: {n_molova_iz_coolpropa}")


# ================ UKUPNI ZBROJEVI ================
# ================================= m3, Pa i K =================================
# V_LOOP: 0.05
# P_init: 110000.00000000001, T_init: 293.15
# P_target: 4000000.0, T_target: 313.15
# ==============================================================================
# n_total molovi u idealnom slucaju pv=nRT: 76.81453300012309
#        Component  Mol [mol]  V_ideal [L]  V_CoolProp [L]     Doses   V_target  Doses adjusted
# 0  CarbonDioxide    69.9012  1548.874300     1539.907300  154.8874  36.729000          3.6729
# 1       Nitrogen     3.0726    68.082400       68.064800    6.8082   2.000000          0.2000
# 2          Argon     0.2304     5.106200        5.102400    0.5106   0.147400          0.0147
# 3          Water     0.1536     0.000003        0.000003    0.0000   0.000003          0.0000
# 4  SulfurDioxide     3.4567    76.592700       75.089500    7.6593   0.166200          0.0166
# Ukupan volumen pri početnim uvjetima (CoolProp): 1688.16 L
# Ukupan cilj. volumen pri 40 bara i 40°C (CoolProp): 39.04 L
# ====================================================

# =========== NORMALIZIRANO NA 0.05 ===========
#        Component  Mol [mol]  V_ideal [L]  V_CoolProp [L]       Doses   V_target  Doses adjusted
# 0  CarbonDioxide  89.519134  1983.569473     1972.085864  198.356909  47.037079        4.703708
# 1       Nitrogen   3.934932    87.189884       87.167344    8.718937   2.561305        0.256130
# 2          Argon   0.295062     6.539267        6.534400    0.653901   0.188768        0.018826
# 3          Water   0.196708     0.000004        0.000004    0.000000   0.000004        0.000000
# 4  SulfurDioxide   4.426831    98.088619       96.163543    9.808900   0.212844        0.021259
# Ukupan volumen pri početnim uvjetima (CoolProp): 2161.95 L
# Ukupan cilj. volumen pri 40 bara i 40°C (CoolProp): 50.00 L
# =====================================================

# ============= PENG-ROBINSON EOS TLAK ==============
# PR EOS tlak za V=0.05 m3 i T=40°C: 39.81 bara
# =====================================================

# ============= PENG-ROBINSON EOS TLAK ALI NORMALIZED za 0.05 m3 ==============
# PR EOS tlak za V=0.05 m3 i T=40°C: 50.98 bara
# =====================================================