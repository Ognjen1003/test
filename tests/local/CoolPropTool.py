import pandas as pd
import os
import json
from CoolProp.CoolProp import PropsSI

# Učitaj CSV (preskoči drugi red s jedinicama)
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
csv_path = os.path.join(base_dir, 'data', 'coolprop_fluid_properties_with_units.csv')

df = pd.read_csv(csv_path, skiprows=[1])

# Ukloni prazne redove ako ih ima
df.dropna(subset=["Fluid"], inplace=True)

# Za svaki fluid, dohvatimo Tcrit i Pcrit iz CoolProp
def get_critical_props(fluid):
    try:
        tcrit = PropsSI("TCRIT", fluid)
        pcrit = PropsSI("PCRIT", fluid)
        return pd.Series({"Tcrit_K": tcrit, "Pcrit_Pa": pcrit})
    except Exception as e:
        print(f"Greška za {fluid}: {e}")
        return pd.Series({"Tcrit_K": None, "Pcrit_Pa": None})

# Dodaj nove kolone
df[["Tcrit_K", "Pcrit_Pa"]] = df["Fluid"].apply(get_critical_props)

# Spremi kao JSON
df.to_json("coolprop_data_with_critical.json", orient="records", indent=2)
print("✔️ Podaci spremljeni u 'coolprop_data_with_critical.json'")