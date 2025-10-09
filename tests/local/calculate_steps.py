from src.Endpoints.FlowModul import calculate_steps
from src.Models.FlowModels import FlowInputModel
from data.testData import ComponentData
from src.Classes import UtilClass
import src.EnumsClasses.MethodsAndTypes as MT
import json
import pprint  


title_primer = "oxyfuel_comp1_original"  # za prikaz vise, nije elementarno
data_source = ComponentData.oxyfuel_comp1_original # podaci koji se actually prikazuju  
UtilClass.check_total_fraction(data_source, title_primer)


# Inicijalizacija modela s podacima
input_data1 = FlowInputModel()
input_data1.nsteps = 3
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
input_data.case = "pvt"
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
    composition=data_source
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
    composition=data_source
)


# Ispis rezultata
print("\n=== Rezultati izračuna ===")
data1 = json.loads(result1)
# ispis po stepovima
for item1 in data1:
    print(f"\nStep {item1['step']}:")
    pprint.pprint(item1)

print("\n=== Rezultati izračuna ===")
data = json.loads(result)
# ispis po stepovima
for item in data:
    print(f"\nStep {item['step']}:")
    pprint.pprint(item)
