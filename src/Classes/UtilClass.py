import numpy as np
from scipy.optimize import fsolve
from numpy.ma.core import log10
import pandas as pd
from scipy.interpolate import griddata
import CoolProp.CoolProp as CP
from numba import jit
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

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
                c.fraction,
                c.CpA,
                c.CpB,
                c.CpC,
                c.CpD
            )
        )
    return cpp_list     

@staticmethod
def display(temperatures, pressures, elapsed_time, results_float, title):
    custom_cmap = sns.color_palette("viridis", as_cmap=True)
    data_masked = results_float.mask((results_float == 0) | (results_float == 1), np.nan)


    cmap_combined = ListedColormap(["black", "white"])
    bounds = [0, 0.5, 1.5]
    norm = BoundaryNorm(bounds, cmap_combined.N)

    plt.figure(figsize=(14, 8))

    ax = sns.heatmap(
    data_masked,
    cmap=custom_cmap,
    cbar=True,
    cbar_kws={'label': '[V]'},
    xticklabels=10,
    yticklabels=False 
)

#   Postavi y tickove ručno na svaki 10 bar (npr. 30, 40, ..., 90)
    pressures_rounded = pressures.round(1)  # ako imaš 30.0, 30.1, ...
    tick_values = list(range(int(pressures[0]), int(pressures[-1])+1, 10))  # 30, 40, ...
    tick_indices = [np.abs(pressures_rounded - val).argmin() for val in tick_values]  # indeksi najbliži traženim vrijednostima
    ax.set_yticks(tick_indices)
    ax.set_yticklabels(tick_values)

    temperatures_rounded = temperatures.round(1)
    x_tick_values = list(range(int(temperatures[0]), int(temperatures[-1])+1, 5))
    x_tick_indices = [np.abs(temperatures_rounded - val).argmin() for val in x_tick_values]
    ax.set_xticks(x_tick_indices)
    ax.set_xticklabels(x_tick_values, rotation=45)

    plt.gcf().patch.set_facecolor('#f0f0f0')
    plt.gca().set_facecolor('#f0f0f0')
    plt.gca().invert_yaxis()

    plt.xlabel("Temperatura [K]")
    plt.ylabel("Tlak [bar]")
    plt.title(title)

    print(f"Vrijeme : {elapsed_time:.5f} sek")

    plt.tight_layout()
    plt.show()

