from src.Classes import UtilClass
from data.testData import ComponentData
import pandas as pd

# učitavanje excela  ==============   case1    ===============

# pvt_verzija = True

# if pvt_verzija:
#     df = pd.read_excel("data\excel\lookup_table_case1_pvt_version.xlsx", engine="openpyxl")
#     results = df.pivot(index="P_bar", columns="T_K", values="Fv")
#     title = "Lookup_table case1 PVT"
# else:
#     df = pd.read_excel("data\excel\lookup_table_case1.xlsx", engine="openpyxl")   
#     results = df.pivot(index="p", columns="t", values="Fv")
#     title = "Lookup_table case1 file"


# results_display = results.astype(float)
# temperatures = results.columns.to_list()
# pressures = results.index.to_list()

# UtilClass.display(temperatures, pressures, results_display, title, False)


# učitavanje excela  ==============   case1    ===============

title = "Zebanec study 2007"

UtilClass.check_total_fraction(ComponentData.data_zebanec2_study_2007, title)