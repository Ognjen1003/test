from Models.Component import Component
from data.testData import ComponentData
from src.Classes.EOS.PREnthalpyCalc import PREnthalpyCalc

#ovo je problem jer kod single phase slucaja RachfordRice.solve vrati -1
components = []

for comp in ComponentData.example_flow_loop:
    component = Component(
        name=comp.name,
        formula=comp.formula,
        Mw=comp.Mw,
        Tc=comp.Tc,
        Pc=comp.Pc,
        omega=comp.omega,
        fraction=comp.fraction,
        CpA=comp.CpA,
        CpB=comp.CpB,
        CpC=comp.CpC,
        CpD=comp.CpD
    )
    components.append(component)


T1=270
P1=45

T1=320
P1=65

enthalpy_calc = PREnthalpyCalc(components)  
total_enthalpy1 = enthalpy_calc.compute_enthalpy_from_z(T1, P1)
total_enthalpy2 = enthalpy_calc.compute_enthalpy_from_z(T1, P1)

print(f"{total_enthalpy1} ---- {total_enthalpy2}")
print(f"===========================")
print(f"{total_enthalpy1 -total_enthalpy2}")
