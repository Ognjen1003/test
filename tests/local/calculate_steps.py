from src.Endpoints.FlowModul import calculate_steps
from src.Models.FlowModels import FlowInputModel
from data.testData import ComponentData
from src.Classes import UtilClass
import src.EnumsClasses.MethodsAndTypes as MT
import json
import pprint  

# Inicijalizacija modela s podacima
input_data = FlowInputModel()
input_data.nsteps = 1
input_data.L = 40000
input_data.d_in = 0.315925
input_data.e = 0.0001
input_data.p = 410000
input_data.T = 275
input_data.qm = 23.75
input_data.case = "case1"
input_data.visual = 0  # ako je dio modela


title_primer = "oxyfuel_comp1"  # za prikaz vise, nije elementarno
data_source = ComponentData.oxyfuel_comp1 # podaci koji se actually prikazuju  

UtilClass.check_total_fraction(data_source, title_primer)


# Poziv funkcije
result = calculate_steps(
    steps=input_data.nsteps,
    length=input_data.L,
    d_in=input_data.d_in,
    e=input_data.e,
    p=input_data.p,
    tK=input_data.T,
    qm=input_data.qm,
    case=input_data.case,
    composition=data_source
)


# Ispis rezultata
print("\n=== Rezultati izračuna ===")
data = json.loads(result)
# ispis po stepovima
for item in data:
    print(f"\nStep {item['step']}:")
    pprint.pprint(item)
