import src.EnumsClasses.MethodsAndTypes as MT
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
from typing import List


bin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
sys.path.insert(0, bin_path)  # stavim ga kao prvi prioritet
import eos_cpp 


class Util:
    def __init__(self):
        pass

    @staticmethod
    def check_total_fraction(data, label):
        total = sum(comp.fraction for comp in data)
        if abs(total - 1.0) < 1e-6:
            print(f"{label}: OK (sum = {total:.6f})")
        else:
            print(f"{label}: NOT OK (sum = {total:.6f})")
            raise ValueError("fraction sum not 1.0")

    @staticmethod
    def convert_to_cpp_components(py_components):
        cpp_list = []
        for c in py_components:
            cpp_list.append(
                eos_cpp.Component(
                    c.name,
                    c.formula,
                    c.Mw,
                    c.Tc,
                    c.Pc,
                    c.omega,
                    c.fraction
                )
            )
        return cpp_list     

    @staticmethod
    def display(temperatures, pressures, results_float, title, display_iterations):
        custom_cmap = sns.color_palette("viridis", as_cmap=True)

        arr = np.asarray(results_float, dtype=float)
        nP, nT = len(pressures), len(temperatures)

        if arr.shape == (nT, nP):
            arr = arr.T  
        elif arr.shape != (nP, nT):
            raise ValueError(f"results_float shape {arr.shape} ne odgovara (len(pressures), len(temperatures)) = {(nP, nT)}")

        # --- maskiraj sentinel vrijednosti ---
        arr = np.where((arr == -1), np.nan, arr)

        # --- okreni redke tako da nizak tlak bude DOLJE (bez invert_yaxis) ---
        arr_plot = arr[::-1, :]                # flip po y
        pressures_plot = np.asarray(pressures)[::-1]

        plt.figure(figsize=(14, 8))
        ax = sns.heatmap(
            arr_plot,
            cmap=custom_cmap,
            cbar=True,
            cbar_kws={'label': '[V]'},
            xticklabels=False,   # ručno ćemo postaviti tickove
            yticklabels=False
        )

        # --- opcionalni overlay iteracija (1..5 u nijansama roze) ---
        if display_iterations:
            pinks = {
                1: "#ffd1dc",
                2: "#ffb3c6",
                3: "#ff8fab",
                4: "#f15bb5",
                5: "#c9184a",
            }
            for val, color in pinks.items():
                overlay = np.where(np.isclose(arr_plot, float(val)), 1.0, np.nan)
                cmap_one = ListedColormap([color])
                cmap_one.set_bad(alpha=0.0)
                ax.imshow(
                    overlay,
                    cmap=cmap_one,
                    interpolation='nearest',
                    origin='upper',   # sad je sve već okrenuto
                    aspect='auto',
                    alpha=0.9,
                    zorder=3
                )
            legend_patches = [Patch(facecolor=c, edgecolor='k', label=f'={v}') for v, c in pinks.items()]
            ax.legend(handles=legend_patches, title='Posebne vrijednosti', loc='upper right', frameon=True)

        pressures_plot_rounded = pressures_plot.round(1)
        y_vals = list(range(int(min(pressures_plot)), int(max(pressures_plot)) + 1, 10))
        y_idx = [int(np.abs(pressures_plot_rounded - val).argmin()) for val in y_vals]
        ax.set_yticks(y_idx)
        ax.set_yticklabels(y_vals)

        temps = np.asarray(temperatures)
        temps_round = temps.round(1)
        x_vals = list(range(int(temps.min()), int(temps.max()) + 1, 5))
        x_idx = [int(np.abs(temps_round - val).argmin()) for val in x_vals]
        ax.set_xticks(x_idx)
        ax.set_xticklabels(x_vals, rotation=45)

        fig = plt.gcf()
        fig.patch.set_facecolor('#f0f0f0')
        ax.set_facecolor('#f0f0f0')

        plt.xlabel("Temperatura [K]")
        plt.ylabel("Tlak [bar]")
        plt.title(title)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def _mw_array_kg_per_mol(components: List) -> np.ndarray:
        #Vrati molekulske mase u kg/mol; ako su dane u g/mol, konvertira
        M = np.array([getattr(c, "Mw") for c in components], dtype=float)
        M = M / 1000.0
        return M

    @staticmethod
    def _normalize(x: np.ndarray) -> np.ndarray:
        x = np.array(x, dtype=float).copy()
        s = x.sum()
        if s <= 0:
            raise ValueError("Sastav faze ima sumu ≤ 0.")
        return x / s

    @staticmethod
    def _name_key(c) -> str:
        nm = getattr(c, "name", "").upper().replace("₂", "2")
        if "CO2" in nm: return "CO2"
        if "N2"  in nm: return "N2"
        raise ValueError(f"Komponenta '{nm}' nije podržana (očekujem CO2/N2).")

    @staticmethod
    def _Kay_mixing(x: np.ndarray, prop: np.ndarray) -> float:
        """Jednostavno Kay miješanje: sum(x_i * prop_i)."""
        return float(np.dot(x, prop))

    @staticmethod
    def _Vc_from_TcPcZc(Tc: np.ndarray, Pc: np.ndarray, Zc: np.ndarray) -> np.ndarray:
        """Kritični molarni volumen iz Tc,Pc,Zc: Vc = Zc R Tc / Pc  [m^3/mol]."""
        return Zc * MT.CONSTANTS.R * Tc / Pc

    def is_near_critical(Zl: float, Zv: float, tol: float = 1e-6) -> bool:

        # Zl manji  (liquid-like)
        # Zv veći  (vapor-like)

        if Zl is None or Zv is None:
            return False
        return abs(Zv - Zl) < tol

    @staticmethod
    def get_phase(phase: float) -> MT.Phase:

        if phase == -2 or phase == -3 or phase == 2 or phase == -10:                        #tekuce
            return MT.Phase.LIQUID
        elif phase == -8 or phase == -9 or phase == 3 or phase == 5:                       #plinovito
            return MT.Phase.VAPOR
        
        return MT.Phase.VAPORLIQUID