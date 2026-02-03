from src.Endpoints.FlowModul import calculate_steps2
from src.Models.FlowModels import FlowInputModel
from data.testData import ComponentData
from src.Classes.UtilClass import Util
from src.Models.Fitting import FittingK
import src.EnumsClasses.MethodsAndTypes as MT
import pandas as pd


title_primer = "oxyfuel_comp1_original"  # za prikaz vise, nije elementarno
data_source = ComponentData.oxyfuel_comp1_original # podaci koji se actually prikazuju  
Util.check_total_fraction(data_source, title_primer)

#svaka funkcija bitna za izracun ima svojeg blizanca

base = FlowInputModel()
base.L = 40000
base.d_in = 0.315925
base.e = 0.0001
base.p = 4100000
base.T = 275
base.qm = 23.75
base.case = "case1" 
base.visual = 0

input_a = base.model_copy(deep=True)
input_a.fittings = []                         # nema koljena
input_a.step_mode = MT.StepMode.BY_FITTINGS      # nije bitno u ovom slučaju
input_a.max_step_m = None                     # ne koristimo refinement
input_a.virtual_steps = 8                     # <-- KLJUČNO: 8 segmenata

""" result_a = calculate_steps2(
    length=input_a.L,
    d_in=input_a.d_in,
    e=input_a.e,
    p=input_a.p,
    tK=input_a.T,
    qm=input_a.qm,
    case=input_a.case,
    fittings=input_a.fittings,
    step_mode=input_a.step_mode,
    max_step_m=input_a.max_step_m,
    virtual_steps=input_a.virtual_steps,
) """


input_b = base.model_copy(deep=True)
input_b.virtual_steps = None                  # nema virtual koraka
input_b.step_mode = MT.StepMode.BY_FITTINGS
input_b.max_step_m = 5000                     # (za sad) bez refinementa
input_b.fittings = [
    FittingK(at_m=10000.0, K=0.9, kind="elbow"),
    FittingK(at_m=17000.0, K=1.2, kind="elbow"),
    FittingK(at_m=19000.0, K=5.0, kind="valve"),
    FittingK(at_m=20000.0, K=1.2, d_in_new=0.250, kind="contraction"),
    FittingK(at_m=34000.0, K=0.3, d_in_new=0.315925, kind="expansion"),
    FittingK(at_m=37000.0, K=0.3, d_in_new=0.325925, kind="expansion", K_ref="DOWNSTREAM")
]

result_b = calculate_steps2(
    length=input_b.L,
    d_in=input_b.d_in,
    e=input_b.e,
    p=input_b.p,
    tK=input_b.T,
    qm=input_b.qm,
    case=input_b.case,
    fittings=input_b.fittings,
    step_mode=input_b.step_mode,
    max_step_m=input_b.max_step_m,
    virtual_steps=input_b.virtual_steps,
)



def print_res(result):
    cols_round = ["p1", "p2", "dp_fric", "dp_minor", "dp_total", "u", "Re", "ff", "rho_g", "mu"]
    df = pd.DataFrame(result)  # result = ona lista dictova

    for c in cols_round:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.round({
    "L": 0,
    "p1": 3,
    "p2": 3,
    "dp_fric": 3,
    "dp_minor": 3,
    "dp_total": 3,
    "u": 3,
    "Re": 0,
    "ff": 6,
    "rho_g": 3,
    "mu": 10,
    "K_step": 3 })

    print(df.to_string(index=False))


print("resultA - nema koljena - 8 koraka")
print_res(result_a)

print("resultB - 2 koljena - 3 koraka")
print_res(result_b)

