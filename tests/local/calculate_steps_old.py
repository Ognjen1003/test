from src.Endpoints.FlowModul import calculate_steps
from src.Models.FlowModels import FlowInputModel
from data.testData import ComponentData
from src.Classes.UtilClass import Util
import src.EnumsClasses.MethodsAndTypes as MT
import pandas as pd


title_primer = "oxyfuel_comp1_original"  # za prikaz vise, nije elementarno
data_source = ComponentData.oxyfuel_comp1_original # podaci koji se actually prikazuju  
Util.check_total_fraction(data_source, title_primer)

def print_res(result):
    cols_round = ["p1", "p2", "dp_fric", "dp_minor", "dp_total", "u", "Re", "ff", "rho_g", "mu"]
    df = pd.DataFrame(result)  # result = ona lista dictova

    for c in cols_round:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.round({
    "step": 0,
    "L": 3,
    "p1": 3,
    "t": 3,
    "mu": 3,
    "rho_g": 3,
    "u": 3,
    "Re": 0,
    "ff": 6,
    "rho_g": 3,
    "dp": 10,
    "p2": 3 })

    print(df.to_string(index=False))

# Inicijalizacija modela s podacima
input_data1 = FlowInputModel()
input_data1.nsteps = 8
input_data1.L = 40000
input_data1.d_in = 0.315925
input_data1.e = 0.0001
input_data1.p = 4100000
input_data1.T = 275
input_data1.qm = 23.75
input_data1.case = "case1"
input_data1.visual = 0  # ako je dio modela

# Inicijalizacija modela s podacima
input_data = FlowInputModel()
input_data.nsteps = 3
input_data.L = 40000
input_data.d_in = 0.315925
input_data.e = 0.0001
input_data.p = 4100000
input_data.T = 275
input_data.qm = 23.75
input_data.case = "oxy1"
input_data.visual = 0  # ako je dio modela





# Poziv funkcije
result1 = calculate_steps(
    steps=input_data1.nsteps,
    length=input_data1.L,
    d_in=input_data1.d_in,
    e=input_data1.e,
    p=input_data1.p,
    tK=input_data1.T,
    qm=input_data1.qm,
    case=input_data1.case,
)

result = calculate_steps(
    steps=input_data.nsteps,
    length=input_data.L,
    d_in=input_data.d_in,
    e=input_data.e,
    p=input_data.p,
    tK=input_data.T,
    qm=input_data.qm,
    case=input_data.case,
)


# Ispis rezultata
print("\n=== Rezultati izračuna ===")
print_res(result1)

print("\n=== Rezultati izračuna ===")
print_res(result)