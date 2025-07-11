from src.Endpoints.FlowModul import calculate_steps
from src.Models.FlowModels import FlowInputModel
import pprint  

# Inicijalizacija modela s podacima
input_data = FlowInputModel()
input_data.nsteps = 2
input_data.L = 40000
input_data.d_in = 0.315925
input_data.e = 0.0001
input_data.p = 4000000
input_data.T = 293.15
input_data.qm = 23.75
input_data.case = "case1"
input_data.visual = 0  # ako je dio modela

# Poziv funkcije
result = calculate_steps(
    steps=input_data.nsteps,
    length=input_data.L,
    d_in=input_data.d_in,
    e=input_data.e,
    p=input_data.p,
    tK=input_data.T,
    qm=input_data.qm,
    case=input_data.case
)

# Ispis rezultata
print("\n=== Rezultati izračuna ===")
for i, row in enumerate(result):
    print(f"\nStep {i+1}:")
    pprint.pprint(row)
