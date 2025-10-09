from typing import List, Dict
import src.Classes.EOS.RachfordRice as RR
from src.Classes.EOS import PengRobinsonEOS
from Models.Component import Component
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import pandas as pd
import src.EnumsClasses.MethodsAndTypes as MT

# ===== GENERALNE STVARI ZA SAD =====
R = 8.314462618  # J/(mol·K)
V_LOOP = 0.05  # m3
P_target_bar=60
T_target_C=60
P_init_bar=1.1
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
    # pretpostavka da je ovo plin uvijek t da znam gdje sam u phase diagramu
    # Phase = RachfordRice.return_phase()  
    #if....
    #else...
    p_real = eos.get_pressure(v_molar)
    return p_real


if __name__ == "__main__":
    
    print("\n================ UKUPNI ZBROJEVI ================")
    df = calculate_charging()
    print(df)

    total_coolprop = df["V_CoolProp [L]"].astype(float).sum()
    total_target = pd.to_numeric(df["V_target"], errors='coerce').sum()

    print(f"Ukupan volumen pri početnim uvjetima (CoolProp): {total_coolprop:.2f} L")
    print(f"Ukupan cilj. volumen pri {P_target_bar} bara i {T_target_C}°C (CoolProp): {total_target:.2f} L")
    print("====================================================")

    
    print(f"\n=========== NORMALIZIRANO NA {V_LOOP} ===========")
    df_normalized = normalize_mols_to_volume(df, target_volume_liters=50)
    print(df_normalized)

    total_coolprop = df_normalized["V_CoolProp [L]"].astype(float).sum()
    total_target = pd.to_numeric(df_normalized["V_target"], errors='coerce').sum()

    print(f"Ukupan volumen pri početnim uvjetima (CoolProp): {total_coolprop:.2f} L")
    print(f"Ukupan cilj. volumen pri {P_target_bar} bara i {T_target_C}°C (CoolProp): {total_target:.2f} L")
    print("=====================================================")

    print("\n============= PENG-ROBINSON EOS TLAK ==============")
    T_K = T_target_C + 273.15
    pressure_real = calculate_pressure_from_eos(df, T_K, V_LOOP)
    print(f"PR EOS tlak za V={V_LOOP} m3 i T={T_target_C}°C: {pressure_real / 1e5:.2f} bara")
    print("=====================================================")

    print(f"\n============= PENG-ROBINSON EOS TLAK ALI NORMALIZED za {V_LOOP} m3 ==============")
    T_K = T_target_C + 273.15
    pressure_real = calculate_pressure_from_eos(df_normalized, T_K, V_LOOP)
    print(f"PR EOS tlak za V={V_LOOP} m3 i T={T_target_C}°C: {pressure_real / 1e5:.2f} bara")
    print("=====================================================")