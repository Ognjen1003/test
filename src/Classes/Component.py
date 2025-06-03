from pydantic import BaseModel

class Component(BaseModel):
    name: str
    formula: str
    Mw: float
    Tc: float
    Pc: float
    omega: float
    fraction: float
    CpA: float
    CpB: float
    CpC: float
    CpD: float
